
library(equivalence)

get_symbols <- function(df, feature){

	df <- as.data.frame(df)
	if (feature == "transcript"){
		rownames(df) <- df[,1]
		df <- df[,-1]
		symbols <- lookup_symbols_from_transcripts(rownames(df))
		rownames(df) <- str_replace(rownames(df), "\\.\\d+.*", "")
		df <- merge(symbols, df, by = 0)
	} else if (feature == "gene"){
		rownames(df) <- df[,1]
		df <- df[,-1]
		symbols <- lookup_symbols_from_genes(rownames(df))
		df <- tibble::rownames_to_column(df, var = "gene_id")
		symbols <- tibble::rownames_to_column(symbols, var = "gene_id")
		df["gene_id"] <- str_replace(df["gene_id"], "\\.\\d+.*", "")
		browser()
		df <- merge(symbols, df, by = "gene_id")
	}
	
	return(df)
}

name_disease_key <- as.character(sample_sheet$disease_group)
names(name_disease_key) <- sample_sheet$study_id	

counts2 <- tidyr::gather(counts, "sample", "counts", 2:25) %>% 
	dplyr::mutate(disease_group = recode(sample, !!!name_disease_key)) %>%
	tidyr::spread(disease_group, counts, fill = 0) %>%
	group_by(gene_id) %>%
	identity()

test <- counts2 %>% 
	summarize(ctrl = sum(ctrl), os = sum(os)) %>% 
	dplyr::arrange(desc(os))

if (feature == "gene"){
	test <- as.data.frame(test)
	
	rownames(test) <- test[,1]
	test <- test[,-1]
	symbols <- lookup_symbols_from_genes(rownames(test))
	symbols <- tibble::rownames_to_column(symbols, var = "gene_id")
	test <- rownames_to_column(test, var = "gene_id")
	test$gene_id <- str_replace(test$gene_id, "\\.\\d+.*", "")
	test <- merge(symbols, test, by = "gene_id")
}

test <- dplyr::arrange(test, desc(os))

counts3 <- counts2 %>% 
	split(.$gene_id) %>% 
	identity()

sum_dfs <- purrr::map(counts3, function(x) colSums(x[c("ctrl", "os")]))

all_zeros <- which(unlist(purrr::map(sum_dfs, function(x) all(x != 0))))

counts3 <- counts3[all_zeros]

test2 <- counts3[1]

tost_res <- purrr::map(counts3, function(x) equivalence::tost(x[["ctrl"]], x[["os"]]))

tidy_tost <- function(tost_out, name){
	feature = name
	mean.x = tost_out$estimate[1]
	mean.y = tost_out$estimate[2]
	df <- tost_out$parameter
	tost.p.value <- tost_out$tost.p.value
	result <- tost_out$result
	
	return(list("feature" = feature, "mean.x" = mean.x, "mean.y" = mean.y, "df" = df, "p.value" = tost.p.value, "result" = result))
	
}

tost_df <- purrr::map2_df(tost_res, names(tost_res), tidy_tost) %>% 
	dplyr::arrange(desc(mean.y)) %>% 
	# dplyr::filter(result == "rejected") %>% 
	identity()

tost_df <- as.data.frame(tost_df)

if (feature == "transcript"){
	rownames(tost_df) <- tost_df[,1]
	tost_df <- tost_df[,-1]
	symbols <- lookup_symbols_from_transcripts(rownames(tost_df))
	rownames(tost_df) <- str_replace(rownames(tost_df), "\\.\\d+.*", "")
	tost_df <- merge(symbols, tost_df, by = 0)
} else if (feature == "gene"){
	rownames(tost_df) <- tost_df[,1]
	tost_df <- tost_df[,-1]
	symbols <- lookup_symbols_from_genes(rownames(tost_df))
	rownames(tost_df) <- str_replace(rownames(tost_df), "\\.\\d+.*", "")
	tost_df <- merge(symbols, tost_df, by = 0)
}

tost_df <- dplyr::arrange(tost_df, desc(mean.y)) %>% 
	# dplyr::filter(result == "rejected") %>% 
	identity()

