library(tidyverse)
library(stringdist)

# Get the keystreams from PDF and prepend 0
keystreams <- c("0002000010110102111112122210011122221010102121222022000221201020221002121121000212222021211121122221",
                "0202020122121210120001200210222112020222022222220220001221012111022121120202022211221112202002121022",
                "0221221101200221120220011002222100000020200021121021020122100021201010210202002000101020022121100100",
                "0100122100011112100120210020011102201122122100100120122212000021220022012202201100010212222110222020")

# Function to run provided perl script
runPerlscript <- function(DNA, keystream){
  command <- str_c("./keystream_randomization.pl ",
                   DNA,
                   " ",
                   keystream,
                   sep = "")
  
  result <- system(command, intern = T)
  
  return(result)
}

seqs <- read_tsv("DNA/IDcheck_out.tsv") %>%
  mutate(mod4 = index %% 4) %>%
  mutate(keystream = keystreams[mod4 + 1]) %>% 
  mutate(decodedDNA = map2_chr(code, keystream, ~ runPerlscript(.x, .y))) %>%
  mutate(decodedDNA = if_else(orientation == "REVERSE", # Reverse complement of reverse orientation
         stringi::stri_reverse(chartr("ACTG", "TGAC", decodedDNA)),
         decodedDNA))

seqs %>%
  select(seqname, decodedDNA) %>%
  write_tsv("DNA/all_decoded_seqs.tsv", col_names = F)

# Plot
seqs %>%
  ggplot(aes(x = ID)) +
  geom_bar()

seqs %>%
  select(ID, index, count, decodedDNA ) %>%
  write_tsv("seqs.tsv")


# save scores
scores <- tibble(positions = as.character(),
                 coverage = as.integer(),
                 sample = as.integer())

# Loop over all samples, except 7
for (j in c(0:6,8)){
  
  IDj <- seqs %>%
    filter(ID == j) 
  
  # Initiate empty vector
  concatseq <- ""
  
  # Loop over all indices
  for (i in 0:max(IDj$index)){
    # Index 0
    if (i == 0){
      
      # Identify the first 25 bp of the DNA-storage string
      df <- IDj %>% 
        filter(index == i) %>% # filter only for reads with index 0
        select(count, decodedDNA) %>%
        mutate(seq_to_concat = str_sub(decodedDNA, 1, 25)) %>% # take the first 25 bp of the decodedDNA
        group_by(seq_to_concat) %>%
        add_count() %>% # Figure out how many times this small subpart of the sequence can be found in all seqs
        ungroup() %>%
        mutate(count = count * n) %>% # Multiply with read count to get total amount
        select(-decodedDNA, -n) %>%
        distinct() %>%
        filter(count == max(count)) 
      
      # Add concatseq
      concatseq <- df$seq_to_concat[1]
      
      # Add count to scores
      scores <- scores %>%
        add_row(positions = "0-25", 
                coverage = df$count,
                sample = j)
      
    # Index 1
    } else if ( i == 1 ) {
      
      # Identify bases 26 - 50
      df <- IDj %>% 
        filter(index %in% c(i-1, i)) %>% # filter only for reads with index 0 and index 1 
        select(index, count, decodedDNA) %>%
        mutate(seq_to_concat = if_else(index == i, 
                                       str_sub(decodedDNA, 1, 25), # take the first 25 bp of the decodedDNA for index i
                                       str_sub(decodedDNA, 26, 50))) %>% # But take bp 26 - 50 for index i -1
        group_by(seq_to_concat) %>%
        add_count() %>% # Figure out how many times this small subpart of the sequence can be found in all seqs
        mutate(count = count * n,
               total_count = sum(count)) %>% # Multiply with read count to get total amount
        ungroup() %>%
        select(-decodedDNA, -n, - index, - count) %>%
        distinct() %>%
        filter(total_count == max(total_count)) 

      
      # Add to previous seq
      concatseq <- str_c(concatseq, df$seq_to_concat[1])
      
      # Add count to scores
      scores <- scores %>%
        add_row(positions = str_c((i*25) + 1, (i+1) * 25, sep = "-"), 
                coverage = df$total_count,
                sample = j)
    
    # Index 2    
    } else if ( i == 2 ){
      
      # Identify bases 51 - 75
      df <- IDj %>% 
        filter(index %in% c(i-2, i-1, i)) %>% # filter only for reads with index 0, 1, 2 
        select(index, count, decodedDNA) %>%
        mutate(seq_to_concat = if_else(index == i, 
                                       str_sub(decodedDNA, 1, 25), # take the first 25 bp of the decodedDNA for index i
                                       if_else(index == i - 1, 
                                               str_sub(decodedDNA, 26, 50), # But take bp 26 - 50 for index i -1
                                               str_sub(decodedDNA, 51, 75)))) %>%  # And take bp 51 - 75 for i -2
        group_by(seq_to_concat) %>%
        add_count() %>% # Figure out how many times this small subpart of the sequence can be found in all seqs
        mutate(count = count * n,
               total_count = sum(count)) %>% # Multiply with read count to get total amount
        ungroup() %>%
        select(-decodedDNA, -n, - index, - count) %>%
        distinct() %>%
        filter(total_count == max(total_count)) 
      
      
      # Add to previous seq
      concatseq <- str_c(concatseq, df$seq_to_concat[1])
      
      # Add count to scores
      scores <- scores %>%
        add_row(positions = str_c((i*25) + 1, (i+1) * 25, sep = "-"), 
                coverage = df$total_count,
                sample = j)
      
    # Index 3    
    } else if ( i == 3 ){
      
      # Identify bases 76 - 100
      df <- IDj %>% 
        filter(index %in% c(i-3, i-2, i-1, i)) %>% # filter only for reads with index 0, 1, 2, 3 
        select(index, count, decodedDNA) %>%
        mutate(seq_to_concat = if_else(index == i, 
                                       str_sub(decodedDNA, 1, 25), # take the first 25 bp of the decodedDNA for index i
                                       if_else(index == i - 1, 
                                               str_sub(decodedDNA, 26, 50), # But take bp 26 - 50 for index i -1
                                               if_else(index == i - 2, 
                                                       str_sub(decodedDNA, 51, 75), # Take bp 51 - 75 for i -2
                                                       str_sub(decodedDNA, 76, 100))))) %>%  # And take bp 76 - 100 for i - 3
        group_by(seq_to_concat) %>%
        add_count() %>% # Figure out how many times this small subpart of the sequence can be found in all seqs
        mutate(count = count * n,
               total_count = sum(count)) %>% # Multiply with read count to get total amount
        ungroup() %>%
        select(-decodedDNA, -n, - index, - count) %>%
        distinct() %>%
        filter(total_count == max(total_count)) 
      
      
      # Add to previous seq
      concatseq <- str_c(concatseq, df$seq_to_concat[1])
      
      # Add count to scores
      scores <- scores %>%
        add_row(positions = str_c((i*25) + 1, (i+1) * 25, sep = "-"), 
                coverage = df$total_count,
                sample = j)
      
    } else if (any(c(i-3, i-2 , i-1, i) %in% IDj$index)) { # This awefull condition checks whether at least one read contains information for the wanted positions
      
      # Identify bases (i*25) + 1 until (i+1) * 25
      df <- IDj %>% 
        filter(index %in% c(i-3, i-2, i-1, i)) %>% # filter only for reads with right indices
        select(index, count, decodedDNA) %>%
        mutate(seq_to_concat = if_else(index == i, 
                                       str_sub(decodedDNA, 1, 25), # take the first 25 bp of the decodedDNA for index i
                                       if_else(index == i - 1, 
                                               str_sub(decodedDNA, 26, 50), # But take bp 26 - 50 for index i -1
                                               if_else(index == i - 2, 
                                                       str_sub(decodedDNA, 51, 75), # Take bp 51 - 75 for i -2
                                                       str_sub(decodedDNA, 76, 100))))) %>%  # And take bp 76 - 100 for i - 3
        group_by(seq_to_concat) %>%
        add_count() %>% # Figure out how many times this small subpart of the sequence can be found in all seqs
        mutate(count = count * n,
               total_count = sum(count)) %>% # Multiply with read count to get total amount
        ungroup() %>%
        select(-decodedDNA, -n, - index, - count) %>%
        distinct() %>%
        filter(total_count == max(total_count)) 
      
      
      # Add to previous seq
      concatseq <- str_c(concatseq, df$seq_to_concat[1])
      
      # Add count to scores
      scores <- scores %>%
        add_row(positions = str_c((i*25) + 1, (i+1) * 25, sep = "-"), 
                coverage = df$total_count,
                sample = j)
      
      
    } else {
      break
    }
  }
  
  path <- str_c("DNA/ID", j , ".dna")
  write_file(concatseq, path)
}

scores %>%
  ggplot(aes(x = positions, y = coverage)) +
  geom_point() +
  facet_wrap(~ sample, nrow = 8) 

write_tsv(scores, "scores.tsv")


