def get_fastq_extention(samples):
   ext_set = set([''.join(pathlib.Path(x).suffixes) for x in samples.filename.to_list()])
   if len(ext_set) > 1:
      return "fasta files must have consistent extensions"
   else:
      return ext_set.pop()