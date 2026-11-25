#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Import individual modules
include {staging}        from './modules/staging'
include {illumination}   from './modules/illumination'
include {registration}   from './modules/registration'
include {dearray}        from './modules/dearray'
include {segmentation}   from './modules/segmentation'
include {quantification} from './modules/quantification'
include {downstream}     from './modules/downstream'
include {viz}            from './modules/viz'
include {background}     from './modules/background'

params.in = null

// Helper functions for finding raw images and precomputed intermediates
def findFiles0(pre, key, pattern) {
    pre[key] ?
        channel.fromPath("${params.in}/$key/$pattern") : channel.empty()
}
def findFiles(pre, key, pattern, ife) {
    pre[key] ?
        channel.fromPath("${params.in}/$key/$pattern").ifEmpty(ife) : channel.empty()
}
def findDirs(pre, key, ife) {
    pre[key] ?
        channel.fromPath("${params.in}/$key/*", type: 'dir').ifEmpty(ife) : channel.empty()
}

// Define the primary mcmicro workflow
workflow {
    // Expecting --in parameter
    if( params.in == null )
        error "Please specify the project directory with --in"

    // Parse MCMICRO parameters (mcp)
    def mcp = mcmicro.Opts.parseParams(
        params, 
        "$projectDir/config/schema.yml",
        "$projectDir/config/defaults.yml"
    )

    // Separate out workflow parameters (wfp) to simplify code
    def wfp = mcp.workflow

    // Identify relevant precomputed intermediates
    // The actual paths to intermediate files are given by
    //   pre.collect{ "${params.in}/$it" }
    def pre = mcmicro.Flow.precomputed(wfp)

    // Identify marker information
    // Note: this must stay the first use of `channel` in the workflow. Otherwise,
    // some older Nextflow releases (e.g., 24.10) fail with a compile error instead
    // of reporting the minimum required version.
    def chMrk = channel.fromPath( "${params.in}/markers.csv", checkIfExists: true )

    // Check that deprecated locations are empty
    channel.fromPath( "${params.in}/illumination_profiles/*" )
        .subscribe{ _f ->
        error "illumination_profiles/ is deprecated; please use illumination/ instead"
    }

    // Some image formats store multiple fields of view in a single file. Other
    // formats store each field separately, typically in .tif files, with a separate
    // index file to tie them together. We will look for the index files from
    // multiple-file formats in a first, separate pass in order to avoid finding the
    // individual .tif files instead. If no multi-file formats are detected, then we
    // look for the single-file formats. Also, for multi-file formats we need to
    // stage the parent directory and not just the index file.
    def (formatType, formatPattern) =
        file("${params.in}/raw/**${wfp['multi-formats']}") ?
        ["multi", wfp['multi-formats']] : ["single", wfp['single-formats']]

    def stagingDirs = findDirs(pre, 'staging', 
        {error "No subdirectories found in staging directory"})
    def staging_in = stagingDirs
        .map{ d -> tuple(
            mcmicro.Util.getSampleName(d, file("${params.in}/staging")),
            mcmicro.Util.getCycleNameFromDir(d, file("${params.in}/staging")),
            formatType == "single" ? d : d.parent
        )}
    // Here we assemble tuples of 1) path to stage for each raw image (might be a
    // directory) and 2) relative path to the main file for each image. Processes
    // must input the first as a path and the second as a val to avoid incorrect or
    // redundant file staging. They must also only use the second (relative) path to
    // construct pathnames for scripts etc. mcmicro.Util.escapePathForShell must be
    // used when interpolating these paths into script strings, as we are bypassing
    // the normal way that paths are passed to channels which handles this escaping
    // automatically.
    def rawFiles = findFiles(pre, 'raw', "**${formatPattern}",
        {error "No images found in ${params.in}/raw"})
    def raw = rawFiles
        .map{ f -> tuple(
            mcmicro.Util.getSampleName(f, file("${params.in}/raw")),
            formatType == "single" ? f : f.parent, 
            f
        )}
        .map{ sampleName, toStage, relPath -> 
            tuple(sampleName, toStage, toStage.parent.relativize(relPath).toString()) }

    // Find precomputed intermediates
    def pre_dfp   = findFiles0(pre, 'illumination', "**-dfp.tif")
        .map{ f -> tuple(mcmicro.Util.getSampleName(f, file("${params.in}/illumination")), f) }
    def pre_ffp   = findFiles0(pre, 'illumination', "**-ffp.tif")
        .map{ f -> tuple(mcmicro.Util.getSampleName(f, file("${params.in}/illumination")), f) }
    def pre_img   = findFiles(pre, 'registration', "*.{ome.tiff,ome.tif,tif,tiff,btf,qptiff}",
        {error "No pre-stitched image in ${params.in}/registration"})
    def pre_bsub  = findFiles(pre, 'background', "*.ome.tif",
        {error "No background subtracted image in ${params.in}/background"})
    def pre_bsubm = findFiles(pre, 'background', "*.csv",
        {error "No background subtracted markers file in ${params.in}/background"})
    def pre_cores = findFiles(pre, 'dearray', "*.tif",
        {error "No TMA cores in ${params.in}/dearray"})
    def pre_masks = findFiles(pre, 'dearray', "masks/*.tif",
        {error "No TMA masks in ${params.in}/dearray/masks"})
    def pre_pmap  = findFiles(pre, 'probability-maps', "*/*-pmap.tif",
        {error "No probability maps found in ${params.in}/probability-maps"})
        .map{ f -> tuple(f.getParent().getName(), f) }
        .filter{ t -> wfp['segmentation'].contains(t[0]) }
    def pre_seg   = findFiles(pre, 'segmentation', "**.tif",
        {error "No segmentation masks in ${params.in}/segmentation"})
        .map{ f -> tuple(f.getParent().getName(), f) }.groupTuple()
    def pre_qty   = findFiles(pre, 'quantification', "*.csv",
        {error "No quantification tables in ${params.in}/quantification"})

    // Write out parameters used
    // Note: `params` and `workflow` are null inside the handler, so everything it
    // needs from them must be captured in local variables.
    def path_qc = "${params.in}/qc"
    def roadieImage = params.roadie
    def wfMetadata = workflow
    workflow.onComplete {
        // Create a provenance directory
        file(path_qc).mkdirs()

        // Write out MCMICRO parameters
        def style = new org.yaml.snakeyaml.DumperOptions();
        style.setPrettyFlow(true);
        style.setDefaultFlowStyle(org.yaml.snakeyaml.DumperOptions.FlowStyle.BLOCK);
        file("${path_qc}/params.yml").withWriter{ out -> 
            new org.yaml.snakeyaml.Yaml(style).dump(mcp, out) 
        }

        // Store additional metadata
        file("${path_qc}/metadata.yml").withWriter{ out ->
            out.println "githubTag: $wfMetadata.revision";
            out.println "githubCommit: $wfMetadata.commitId";
            out.println "roadie: $roadieImage";
        }
    }

    staging(mcp, staging_in, chMrk)
    //staging.out.view()
    staging.out.map{
        sample, cycle, path ->
        tuple(sample, cycle, path, path.toString().split('/').last())
        }.toSortedList { a, b -> a[1] <=> b[1] }
        .flatMap()
        .map{
            sample, _cycle, path, name ->
            tuple(sample, path, name)
        }.set{ sorted_staging }
    raw = raw.mix(sorted_staging)

    illumination(wfp, mcp.modules['illumination'], raw)
    registration(mcp, raw,
		 illumination.out.ffp.mix( pre_ffp ),
		 illumination.out.dfp.mix( pre_dfp ))
    img = registration.out.mix(pre_img)

    // Should background subtraction be applied?
    img = img
        .branch{
            nobs: !wfp.background
            bs: wfp.background
        }
    chMrk = chMrk
        .branch{
            nobs: !wfp.background
            bs: wfp.background
        }
    // Apply background if specified
    background(mcp, img.bs, chMrk.bs)
    // Merge against precomputed intermediates
    bsub_image = background.out.image.mix(pre_bsub)
    bsub_marker = background.out.marker.mix(pre_bsubm)
    // Reconcile non-background subtracted and background 
    // subtracted images for downstream processing
    img = img.nobs.mix(bsub_image)
    // Reconcile the marker file to the background subtracted csv
    chMrk = chMrk.nobs.mix(bsub_marker)

    // Are we working with a TMA or a whole-slide image?
    img = img
        .branch {
            wsi: !wfp.tma
            tma: wfp.tma
        }

    // Apply dearray to TMAs only
    dearray(mcp, img.tma)

    // Merge against precomputed intermediates
    tmacores = dearray.out.cores.mix(pre_cores)
    tmamasks = dearray.out.masks.mix(pre_masks)

    // Reconcile WSI and TMA processing for downstream segmentation
    allimg = img.wsi.mix(tmacores)
    segmentation(mcp, allimg, tmamasks, pre_pmap)

    // Merge segmentation masks against precomputed ones and append markers.csv
    segMsk = segmentation.out.mix(pre_seg)
    quantification(mcp, allimg, segMsk, chMrk)

    // Spatial feature tables -> cell state calling
    sft = quantification.out.mix(pre_qty)
    downstream(mcp, sft)

    // Vizualization
    viz(mcp, allimg, chMrk)
}
