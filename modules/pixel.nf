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

    // Isolate cell masks and determine their IDs
    id_msks = segmasks.map{ id, msk -> 
        x = id.split('-',2);
        cellmsk = msk.find{ 
            it.getSimpleName().toString().indexOf('cell') > -1 
        }
        tuple(x[1], x[0], cellmsk) 
    }

    // Isolate spatial feature tables corresponding to the cell masks
    sfts.filter{
        it.getName().toString().indexOf('_cell.csv') > -1
    }

    emit:
        sfts
}