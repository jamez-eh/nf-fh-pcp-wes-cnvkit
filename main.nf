nextflow.preview.dsl=2


process access {
        container 'etal/cnvkit'
        label 'small'

        input:
        path reference

        output:
	path("access.${reference}.bed")

        """
        cnvkit.py access ${reference} -o access.${reference}.bed
        """
}


process target {
	container 'etal/cnvkit'
	label 'small'

	input:
	tuple val(kitID), path(bed_file)
	path refFlat

	output:
	tuple val("${kitID}"), path("${kitID}_targets.bed"), path("${kitID}_antitargets.bed")

	"""
	cnvkit.py target ${bed_file} --annotate ${refFlat} --split --avg-size 200 -o ${kitID}_targets.bed
	cnvkit.py antitarget ${kitID}_targets.bed -o ${kitID}_antitargets.bed

	"""

}


process coverage {
        container 'etal/cnvkit'
        label 'small'

        input:
	tuple val(kitID), path(target), path(antitarget), val(sampleID), val(type), val(patientID), path(bam), val(ploidy), val(purity)
	path access
	path refFlat

        output:
	tuple val("${kitID}"), val("${sampleID}"), val("${type}"), val("${patientID}"), path("${sampleID}.targetcoverage.cnn"), path("${sampleID}.antitargetcoverage.cnn"), val("${ploidy}"), val("${purity}")

        """
        cnvkit.py coverage ${bam} ${target} -o ${sampleID}.targetcoverage.cnn
        cnvkit.py coverage ${bam} ${antitarget} -o ${sampleID}.antitargetcoverage.cnn
        """
}



process pon_reference {
	container 'etal/cnvkit'
        label 'small'

	input:
        tuple val(kitID), path(targets), path(antitargets)
	path(reference)

        output:
	tuple val("${kitID}"), path("${kitID}_reference.cnn")

        """
        echo *
        echo *targetcoverage.cnn
	cnvkit.py reference *targetcoverage.cnn --fasta ${reference} -o ${kitID}_reference.cnn

        """
}


process autobin {
        container 'etal/cnvkit'
        label 'small'
	
	input:
	output:

	"""
	cnvkit.py autobin *.bam -t baits.bed -g access.hg19.bed [--annotate refFlat.txt --short-names]
	"""
}

process	fix {
        container 'etal/cnvkit'
        label 'small'

	input:
	tuple val(kitID), path(pon), val(sampleID), val(type), val(patientID), path(targets), path(antitargets), val(purity), val(ploidy)

	output:
	tuple val("${sampleID}"), val("${patientID}"), path("${sampleID}.cnr"),  val("${purity}"), val("${ploidy}")

	"""
	cnvkit.py fix ${targets} ${antitargets} ${pon} -o ${sampleID}.cnr

	"""
}

process segment {
        container 'etal/cnvkit'
        label 'small'
	
	input:
	tuple val(sampleID), val(patientID), path(cnr), val(purity), val(ploidy)

	output:
	tuple val("${sampleID}"), val("${patientID}"), path("${sampleID}.cns"), val("${purity}"), val("${ploidy}")

	"""
	cnvkit.py segment ${cnr} -o ${sampleID}.cns
	"""
}

process center {
        container 'jamezeh/cnvkit'
        label 'small'

        input:
        tuple val(sampleID), val(patientID), path(cns), val(purity), val(ploidy)

        output:
        tuple val("${sampleID}"), val("${patientID}"), path("${sampleID}_centered.cns"), val("${purity}"), val("${ploidy}")

	"""
	cnvkit.py call -h
	cnvkit.py call -m none ${cns} --center mode -o ${sampleID}_centered.cns

        """

}

process callCNV {
        container 'jamezeh/cnvkit'
        label 'small'

        input:
        tuple val(sampleID), val(patientID), path(cns), val(purity), val(ploidy)

        output:
        tuple val("${sampleID}"), val("${patientID}"), path("${sampleID}_calls.cns")

        """
#	cp -L {vcf} ${sampleID}.vcf.gz
	#gunzip ${sampleID}.vcf.gz 
        cnvkit.py call -t=-1.5,-0.4,0.3,0.7  ${cns} -y -o ${sampleID}_calls.cns 
	#-v {sampleID}.vcf -z 0.25

        """

}

process callCNV_vcf {
        container 'jamezeh/cnvkit'
        label 'small'

        input:
        tuple val(sampleID), val(patientID), path(cns), val(purity), val(ploidy)

        output:
        tuple val("${sampleID}"), val("${patientID}"), path("${sampleID}_calls.cns")

        """
#       cp -L {vcf} ${sampleID}.vcf.gz
        #gunzip ${sampleID}.vcf.gz
        cnvkit.py call -t=-1.5,-0.4,0.3,0.7  ${cns} -y -o ${sampleID}_calls.cns
        #-v {sampleID}.vcf -z 0.25

        """

}





process images {
        container 'etal/cnvkit:0.9.3'
        label 'small'
	
	input:
	tuple val(sampleID), val(patientID), path(cnr), path(cns)
	path(coords)

	output:
        tuple val("${sampleID}"), val("${patientID}"), path("${sampleID}-diagram.pdf"), path("${sampleID}-scatter-genes.pdf"), path("${sampleID}-scatter.pdf")

	
	"""
	cnvkit.py diagram ${cnr} -s ${cns} -o ${sampleID}-diagram.pdf
	cnvkit.py scatter -s *.cn{s,r} -l ${coords} -w 1000000 -o ${sampleID}-scatter-genes.pdf 
	cnvkit.py scatter -s *.cn{s,r} -o ${sampleID}-scatter.pdf

	"""
}


workflow CNVkit_wf {

         take:
         bams_ch
         beds_ch
         reference
	 refFlat
	 
	 main:
	 access(reference)

 	 target(beds_ch, refFlat)

	 bams_target = target.out.cross(bams_ch.map { r-> [r[1], r[0], r[2], r[3], r[4], r[5], r[6]] })
	 	       				    .map{ r -> [r[0][0], r[0][1], r[0][2], r[1][1], r[1][2], r[1][3], r[1][4], r[1][5], r[1][6]]}


	 
	 coverage(bams_target, access.out, refFlat)
	 
	 
  	 normal_coverage = coverage.out.filter { it[2] == 'Normal'}.groupTuple().map { r -> [r[0],r[4], r[5]] }
	 tumor_coverage = coverage.out.filter { it[2] == 'Tumor' | it[2] == 'PDX' }
	 
	 pon_reference(normal_coverage, reference)
	 
	 tumor_coverage.view()
	 pon_crossed = pon_reference.out.cross(tumor_coverage).map{ r -> [r[0][0], r[0][1], r[1][1], r[1][2], r[1][3], r[1][4], r[1][5], r[1][6], r[1][7]] }
	 
	 
	 fix(pon_crossed)
	 fix.out.view()
	 segment(fix.out)
	 center(segment.out)
	
	 callCNV(center.out)
//	 images(fix.out.join(callCNV.out, by : [0,1]), coords)

	 
	 emit:
	 fixed = fix.out 
	 segments = segment.out
	 centered = center.out
	 calls = callCNV.out


}




workflow {
	        
	 beds_ch = Channel
            .fromPath(params.input_beds)
            .splitCsv(header:true)
            .map{ row-> tuple(row.kitID, file(row.capture_bed)) }

	 bams_ch = Channel
	    .fromPath(params.input_csv)
	    .splitCsv(header:true)
	    .map{ row-> tuple(row.sampleID, row.kitID, row.type, row.patientID, file(row.bam), row.ploidy, row.purity) }


	reference = file(params.reference)
	refFlat = file(params.refFlat)
	//coords = file(params.coords)

	 main:

	 CNVkit_wf(bams_ch, beds_ch, reference, refFlat)
	 
	 publish: 
         CNVkit_wf.out.fixed to : "${params.output_folder}/fixed/"
         CNVkit_wf.out.segments to :"${params.output_folder}/segments/"
         CNVkit_wf.out.centered to : "${params.output_folder}/centered/"
         CNVkit_wf.out.calls to : "${params.output_folder}/calls/"


}
