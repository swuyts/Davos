library(tidyverse)
library(oro.dicom)

# Create a function for decoding
fromBasetoTrit <- function(string){
  
  trit <- ""
  
  for (i in 2:nchar(string)){
    
    prevbase <- str_sub(string, i-1, i-1)
    
    triti <- if_else(prevbase == "A",
                     chartr("CGT", "012", str_sub(string,i,i)),
                     if_else(prevbase == "C",
                             chartr("GTA", "012", str_sub(string,i,i)),
                             if_else(prevbase == "G",
                                     chartr("TAC", "012", str_sub(string,i,i)),
                                     if_else(prevbase == "T",
                                             chartr("ACG", "012", str_sub(string, i, i)), 
                                             "FALSE"))))
    trit <- str_c(trit, triti)
    
    
  }
  return(trit)
}

calculate_parity_check <- function(IX_trit){
  str_sub(IX_trit, # extract the necessary elements for the check 
          start = c(1, 3, 5, 7, 9, 11, 13), 
          end = c(1, 3, 5, 7, 9, 11, 13)) %>%
    as.integer() %>% # convert string to number
    sum() %% 3 %>% # calculate sum and modulus
    as.integer()
}

# Read in big df
seqs <- read_tsv("DNA/merged.fastq.assembled.fastq.gz_filtered_trimmed_derep.tsv", col_names = c("seqname", "seq"))

# Split sequence and Identifier 'IX'
seqs <- seqs %>%
  mutate(code = str_sub(seq, 1, 100),
         IX = str_sub(seq, -15, -1)) %>%
  select(-seq)

# Perform QC and figure out ID of reads
IX_check <- seqs %>% 
  mutate(lastnucCode_IX = str_c(str_sub(code, -1, -1), IX)) %>% # Add last nucleotide of seq to IX for converting to trit
  select(-code) %>% # reduce dataset by removing seqs
  mutate(IX_trit = map_chr(lastnucCode_IX, ~fromBasetoTrit(.))) %>% # convert from base to trit
  select(-lastnucCode_IX) %>%
  filter(!str_detect(IX_trit, "A")) %>% # remove seqs where homopolymers are found
  filter(!str_detect(IX_trit, "C")) %>% # remove seqs where homopolymers are found
  filter(!str_detect(IX_trit, "G")) %>% # remove seqs where homopolymers are found
  filter(!str_detect(IX_trit, "T")) %>% # remove seqs where homopolymers are found
  mutate(ID_trit = str_sub(IX_trit, 1 ,2), # Split IX into the three parts
         i3_trit = str_sub(IX_trit, 3, -2),
         P_trit = str_sub(IX_trit,-1,-1))%>%
  mutate(parity_check = map_chr(IX_trit, # Perform parity check 
                            ~calculate_parity_check(.))) %>%
  filter(P_trit == parity_check) %>% # Filter out reads that fail the check
  mutate(ID = strtoi(ID_trit, base = 3), # Convert ID to decimal
         index = strtoi(i3_trit, base = 3), # Convert i3 to decimal
         orientation = if_else(index %% 2 == 0, # Figure out orientation
                              "FORWARD",
                              "REVERSE")) %>%
  select(seqname, ID, index, orientation)

# Filter reads
final <- IX_check %>%
  left_join(seqs) %>%
  mutate(separator = seqname) %>%
  separate(separator, "-", into = c("name", "count")) %>%
  group_by(ID, index) %>%
  arrange(index, ID) 


write_tsv(final, "DNA/IDcheck_out.tsv")

