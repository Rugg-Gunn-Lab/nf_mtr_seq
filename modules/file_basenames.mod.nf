#!/usr/bin/env nextflow
/** This code is from files.mod.nf in the Babraham bioinf nextflow modules. 
It has just been condensed so doesn't search for trimmed files etc.
*/

def getFileBaseNames(fileList) {

    baseNames = [:]

    bareFiles = []

    for (String s : fileList) {
        matcher = s =~ /^(.*)_(R?[1234]).fq.gz$/

        //println (matcher[0])
        if (matcher.matches()) {
            if (! baseNames.containsKey(matcher[0][1])) {
                baseNames[matcher[0][1]] = []
            }
            baseNames[matcher[0][1]].add(matcher[0][2])
        }
        else {
            matcher = s =~ /^(.*).fq.gz$/

            if (matcher.matches()) {
                bareFiles.add(matcher[0][1])
            }
        }
    }
    patterns = []
    for (s in baseNames) {
        //println (s)
        pattern = s.key+"_{"+s.value.join(",")+"}.fq.gz"
        patterns.add(pattern)
        //println("$pattern")
    }
    for (s in bareFiles) {
        pattern = s+".fq.gz"
        patterns.add(pattern)
    }

    return(patterns)
}

// keep this in in case we add renaming to allow fastq.gz files
// def getFileBaseNames(fileList) {

//     baseNames = [:]

//     bareFiles = []

//     for (String s : fileList) {
//         matcher = s =~ /^(.*)_(R?[1234]).(fastq|fq).gz$/

//         //println (matcher[0])
//         if (matcher.matches()) {
//             if (! baseNames.containsKey(matcher[0][1])) {
//                 baseNames[matcher[0][1]] = []
//             }
//             baseNames[matcher[0][1]].add(matcher[0][2])
//         }
//         else {
//             matcher = s =~ /^(.*).(fastq|fq).gz$/

//             if (matcher.matches()) {
//                 bareFiles.add(matcher[0][1])
//             }
//         }
//     }
//     patterns = []
//     for (s in baseNames) {
//         //println (s)
//         pattern = s.key+"_{"+s.value.join(",")+"}.{fastq,fq}.gz"
//         patterns.add(pattern)
//         //println("$pattern")
//     }
//     for (s in bareFiles) {
//         pattern = s+".{fastq,fq}.gz"
//         patterns.add(pattern)
//     }

//     return(patterns)
// }