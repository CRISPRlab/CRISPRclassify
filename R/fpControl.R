filterFalsePositives <- function(df, identity_limit, seq_type) {
  print('...filtering false positives...')
  lociIds <- df %>% distinct(locus)

  exclusionList <- c()
  for (x in 1:nrow(lociIds)){

    if (seq_type == "spacer"){
      ## length test ##
      filterdf <- df %>%
        select(id, sequence, locus, length) %>%
        filter(locus == x)
      avgLength <- filterdf %>%
        group_by(locus) %>%
        summarize(m = mean(length))
      lengthdf <- filterdf %>%
        mutate(percent_shorter = length / avgLength$m) %>%
        filter(percent_shorter < .7)
      if (nrow(lengthdf) > 0){
        exclusionList <- c(exclusionList,avgLength$locus)
        next
      }

      ## identity test ##
      filterdf <- filterdf %>%
        select(id, sequence, locus)
    } else if (seq_type == "repeat"){
      filterdf <- df %>%
        select(id, sequence, locus) %>%
        filter(locus == x)
    } else {
      stop("Invalid seq_type argument passed into filterFalsePositives")
    }

    duplicatedf <- filterdf

    # create cartesian
    cartesiandf <- filterdf %>% inner_join(duplicatedf, by = c("locus")) %>%
      filter(id.x != id.y) # make sure not to measure dist between a repeat and itself

    #calculate distance between each repeat
    cartesiandf <- cartesiandf %>% mutate(identity = percentIden(sequence.x, sequence.y))

    # calc avg identity across entire locus
    avgIden <- cartesiandf %>%
      group_by(locus) %>%
      summarize(m = mean(identity))

    if (seq_type == "repeat"){
      if (avgIden$m < identity_limit) { #.8
        exclusionList <- c(exclusionList,avgIden$locus)
      }
    } else if (seq_type == "spacer"){
      if (avgIden$m > identity_limit) { #.6
        exclusionList <- c(exclusionList,avgIden$locus)
      }
    } else {
      stop("Invalid seq_type argument passed into filterFalsePositives")
    }

  }

  return(df %>% filter(!locus %in% exclusionList))
}
