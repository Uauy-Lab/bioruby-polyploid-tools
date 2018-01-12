module Bio::DB::Blast

	def self.align(opts={})
		target=opts[:target]
		query=opts[:query]
		cmdline = "blastn -query #{query} -db #{target} -outfmt '6 qseqid qstart qend qframe sseqid sstart send sframe score pident	qlen	slen	qseq sseq'"

		status, stdout, stderr = systemu cmdline
		if status.exitstatus == 0
			alns = Array.new unless block_given?
			stdout.each_line do |line|
				puts line
				aln = Bio::DB::ExonerateAlignment.parse_custom(line) 
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


