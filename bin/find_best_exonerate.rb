#!/usr/bin/env ruby


found_cointigs = Set.new
Bio::DB::Exonerate.align({:query=>temp_fasta_query, :target=>target, :model=>model, :chunk=>chunk, :total_chunks=>}) do |aln|
  if aln.identity > min_identity
    exo_f.puts aln.line
    unless found_cointigs.include?(aln.target_id) #We only add once each contig. Should reduce the size of the output file. 
      found_cointigs.add(aln.target_id)
      entry = fasta_file.index.region_for_entry(aln.target_id)
      raise ExonerateException.new,  "Entry not found! #{aln.target_id}. Make sure that the #{target_id}.fai was generated properly." if entry == nil
      region = entry.get_full_region
      seq = fasta_file.fetch_sequence(region)
      contigs_f.puts(">#{aln.target_id}\n#{seq}")
    end
  end  
end