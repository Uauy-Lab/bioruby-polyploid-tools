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
		if fields[7] =="-1"
			fields[7] = "-"
			fields[6] = fields[6].to_i - 1
		else
			fields[7] = "+"
			fields[5] = fields[5].to_i - 1
		end
		fields.join(" ")
	end

	def self.to_vulgar(line)
		qseq, sseq = line.split("\t")[12..13]

		len = qseq.length
		l_status = ""
		l_len = 0
		str = Array.new
		statuses = ""
		for i in 0..len
			if qseq[i] == "-"
				status = "D"
			elsif sseq[i] == "-"
				status = "I"
			else
				status = "M"
			end
			statuses << status
		end
		statuses.split('').each do |e| 
			if l_status != e
				case l_status
				when "M"
					str << ["M", l_len, l_len]
				when "I"
					str << ["G", l_len, 0]
				when "D"
					str << ["G", 0, l_len]
				end
				l_len = 0
			end
			l_status = e
			l_len += 1
		end
		l_len -= 1
		case l_status
		when "M"
			str << ["M", l_len, l_len]
		when "I"
			str << ["G", l_len, 0]
		when "D"
			str << ["G", 0, l_len]
		end

		str.flatten!.join(" ")
	end

	def self.to_exo(line)
		arr = Array.new
		arr << "RESULT:"
		arr << to_sugar(line)
		arr << line.split("\t")[9..11]
		arr << "."
		arr << to_vulgar(line)
		arr.join("\t")
	end

	def self.align(opts={})
		target=opts[:target]
		query=opts[:query]
		max_target_seqs = 15
		max_target_seqs = opts[:max_hits] * 2 if opts[:max_hits]
		cmdline = "blastn -max_target_seqs #{max_target_seqs} -query #{query} -db #{target} -outfmt '6 qseqid qstart qend qframe sseqid sstart send sframe score pident qlen slen qseq sseq'"

		status, stdout, stderr = systemu cmdline
		if status.exitstatus == 0
			alns = Array.new unless block_given?
			stdout.each_line do |e_l|
				#puts e_l
				line = to_exo(e_l)
				#puts line
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
			raise BlasteException.new(), "Error running exonerate. Command line was '#{cmdline}'\n Blast STDERR was:\n#{stderr}"
		end
	end	

	class BlasteException < RuntimeError 
	end	

end


