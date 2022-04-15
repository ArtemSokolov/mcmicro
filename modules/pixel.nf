process cellcutter {
    container "${params.contPfx}${module.container}:${module.version}"

    // Output
    publishDir "${params.in}/cells/${tag}", mode: 'copy', pattern: '*.zip'

    // Provenance
    publishDir "${params.path_prov}", mode: 'copy', pattern: '.command.sh',
        saveAs: {fn -> "${task.name}.sh"}
    publishDir "${params.path_prov}", mode: 'copy', pattern: '.command.log',
        saveAs: {fn -> "${task.name}.log"}

    input:
        val module
        tuple val(tag), path(img), path(msk), path(sft)

    output:
        tuple val(tag), path("*.zip"), emit: cells

    """
    cut_cells -z ${params.cutterOpts} $img $msk $sft .
    """
}

include { getFileID; getImageID } from './lib/util'

workflow pixel {
    take:
        moduleCutter
        imgs
        segmasks
        sfts
    
    main:

    // Determine IDs of images
    id_imgs = imgs.map{ f -> tuple(getImageID(f), f) }

    // Isolate spatial feature tables (SFTs) corresponding to the cell masks
    //  and extract the corresponding IDs
    id_sfts = sfts.filter{
        it.getName().toString().indexOf('_cell.csv') > -1
    }.map{
        f -> tuple(getFileID(f,'_cell'), f)
    }

    // Isolate cell masks and match them up against SFTs and images
    inputs = segmasks.map{ id, msk ->
        cellmsk = msk.find{ 
            it.getSimpleName().toString().indexOf('cell') > -1 
        }
        tuple(id, cellmsk)
    }.combine(id_sfts, by:0)
        .map{ id, msk, sft ->
        x = id.split('-',2); tuple(x[1], id, msk, sft)
    }.combine(id_imgs, by:0)
        .map{ id, tag, msk, sft, img -> tuple(tag, img, msk, sft)}

    // Cut the cells based on combined inputs
    cellcutter(moduleCutter, inputs)

    emit:
        cellcutter.out.cells
}
