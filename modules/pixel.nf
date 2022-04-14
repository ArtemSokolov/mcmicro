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

    // Determine IDs of segmentation masks
    id_msks = segmasks.map{ id, msk -> 
        x = id.split('-',2);
        mskf = msk.find{ it.getSimpleName().toString() == 'cell' }
        tuple(x[1], x[0], mskf) }

    emit:
        id_msks
}