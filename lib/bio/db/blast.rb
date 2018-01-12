module Bio::DB::Blast

	def self.to_sugar(line)
		fields = line.split("\t")[0..8]
		

		if  fields[3] =="-1"
			fields[3] = "-"
			fields[2] = fields[2].to_i - 1
		else
			fields[3] = "+"
			fields[1] = fields[1].to_i - 1
		end

		if  fields[7] =="-1"
			fields[7] = "-"
			fields[6] = fields[6].to_i - 1
		else
			fields[7] = "+"
			fields[5] = fields[5].to_i - 1
		end

		fields.join(" ")
	end

	def self.align(opts={})
		target=opts[:target]
		query=opts[:query]
		cmdline = "blastn -query #{query} -db #{target} -outfmt '6 qseqid qstart qend qframe sseqid sstart send sframe score pident	qlen	slen	qseq sseq'"

		status, stdout, stderr = systemu cmdline
		if status.exitstatus == 0
			alns = Array.new unless block_given?
			stdout.each_line do |line|
				puts line
				arr = line.split("\t")



				aln = Bio::DB::Exonerate::Alignment.parse_custom(line) 
				if aln
					if block_given?
						yield aln
					else
						alns << aln
					end
				end
			end
			return alns unless block_given?
		else
			raise ExonerateException.new(), "Error running exonerate. Command line was '#{cmdline}'\nExonerate STDERR was:\n#{stderr}"
		end
	end	

	class BlasteException < RuntimeError 
	end	

end


