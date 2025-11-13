rm(list = ls())
#Call all used libraries from the beginning
library(tidyverse)
conflicted::conflict_prefer("filter", "dplyr")
library(viridis)
# + scale_color/fill_viridis(discrete = T/F)
theme_set(theme_light())
library("vegan")
library(ggplot2)
library(stringr)
library(dplyr)
library(gplots)
library(VennDiagram)
library(grid)
library(Biostrings)
library(reshape2)
library(ggplot2)
library(foreach)
library(doParallel)

#Read in Onychophora data
dfWORM <- read_tsv("../data/Onychophora.tsv")
problems(dfWORM)
#take a look at the structure of the dataframe
names(dfWORM)

dfWORM %>% 
  filter(!is.na(habitat)) #none of them recorded habitat
dfWORM %>% 
  filter(!is.na(ecoregion)) %>% 
  select(ecoregion)#use ecoregion to get habitats
dfWORM %>% 
  filter(!is.na(ecoregion)) %>% 
  filter(!is.na(species)) %>% 
  select(species) %>% 
  unique() #4 different velvet worm species in this dataset

#Read in Chilophoda data
ddCENT <- read_tsv('../data/Chilopoda.tsv')
#take a look at the structure of the dataframe
names(ddCENT)

ddCENT %>% 
  filter(!is.na(ecoregion)) %>% 
  filter(!is.na(species)) %>% 
  filter(!is.na(biome)) %>% 
  select(species) %>% 
  unique()
ddCENT %>% 
  filter(!is.na(ecoregion)) %>% 
  filter(!is.na(species)) %>% 
  filter(!is.na(biome)) %>% 
  select(biome) %>% 
  unique()

####Create a bar graph showing unique species counts per shared ecoregion
#Number of species in shared ecoregions
VelvetEco <- dfWORM %>% 
  filter(!is.na(ecoregion)) %>% 
  select(ecoregion) %>% 
  unique() %>% 
  print()
CentiEco <- ddCENT %>% 
  filter(!is.na(ecoregion)) %>%
  select(ecoregion) %>% 
  unique() %>% 
  print() 
conflicted::conflicts_prefer(base::intersect)
SharedEco <- intersect(dfWORM$ecoregion, ddCENT$ecoregion) %>% 
  print() #TODO: find out how to remove NA from this list
SharedEco <- SharedEco[!is.na(SharedEco)]
print(SharedEco)

#Step 1: Filter and summarize each dataset
worm_summary <- dfWORM %>%
  filter(ecoregion %in% SharedEco) %>%
  group_by(ecoregion) %>%
  summarise(species_count = n_distinct(species)) %>%
  mutate(group = "Velvet Worm") %>% 
  print()

cent_summary <- ddCENT %>%
  filter(ecoregion %in% SharedEco) %>%
  group_by(ecoregion) %>%
  summarise(species_count = n_distinct(species)) %>%
  mutate(group = "Centipede") %>% 
  print()

# Step 2: Combine into one dataframe
combined_data <- bind_rows(worm_summary, cent_summary)
combined_data <- combined_data %>% 
  mutate(ecoregion = gsub('_', ' ', ecoregion))
print(combined_data)

# Step 3: Plot
ggplot(combined_data, aes(x = str_wrap(ecoregion, width = 15), y = species_count, fill = group)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Unique Species per Shared Ecoregion",
       x = "Ecoregion",
       y = "Number of Unique Species") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        #originally tried angle = 45, hjust = 1, but chao phraya freshwater swamp forests was too long and went right off the page.
        panel.grid.major = element_blank(), #remove major gridlines
        panel.grid.minor = element_blank(),  #remove minor gridlines
        plot.title = element_text(hjust = 0.5)  #center title
  )


####Make a Venn Diagram of total species overlap by ecoregion
#Step 1: isolate ecoregions and species counts for each region
worm_ecoregions <- dfWORM %>%
  group_by(ecoregion) %>%
  summarise(species_count = n_distinct(species)) %>% 
  print()

cent_ecoregions <- ddCENT %>%
  group_by(ecoregion) %>%
  summarise(species_count = n_distinct(species)) %>% 
  print()

#Step 2: Extract unique ecoregions
worm_regions <- worm_ecoregions$ecoregion
worm_regions <- worm_regions[!is.na(worm_regions)] #remove those pesky NAs
print(worm_regions)
cent_regions <- cent_ecoregions$ecoregion
cent_regions <- cent_regions[!is.na(cent_regions)]
print(cent_regions)

#Step 3: clear plots because the function venn.plot overlaps with other graphs
grid.newpage()

#Step 4: Draw 2-set Venn diagram
venn.plot <- draw.pairwise.venn(
  area1 = length(worm_regions),
  area2 = length(cent_regions),
  cross.area = length(intersect(worm_regions, cent_regions)),
  category = c("Onychophora", "Chilopoda"),
  fill = c("purple", "orange"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.5,
  cat.pos = c(0, 180),
  scaled = F
)

grid.text(
  "Ecoregion Overlap by Species", # title text
  y = unit(0.95, "npc"),         # position near top (npc = normalized parent coordinates)
  gp = gpar(fontsize = 16, fontface = "bold")  # text styling
)

####Try and make a pairwise alignment of Onychophora and Chilopoda BIN sequences
#Step 1: Remove NA values, make sure all letters are uppercase to not miss any alignments due to a technicality
dfWORM_clean <- dfWORM %>%
  filter(!is.na(nuc)) %>%
  mutate(nuc = toupper(trimws(nuc)))

ddCENT_clean <- ddCENT %>%
  filter(!is.na(nuc)) %>%
  mutate(nuc = toupper(trimws(nuc)))

#Step 2: Random sampling, because there are too many sequences to do on my computer
#Set a seed for random sampling, for reproducibility purposes
set.seed(123)
#Sample 500 times for both Onychophora and Chilopoda
#500 samples, on my 3 core lenovo thinkpad, took ~45 minutes to complete a pairwise alignment. Feel free to lower the sample amount when running the code to make sure it works.
worm_seq <- sample(DNAStringSet(dfWORM_clean$nuc), 500)
cent_seq <- sample(DNAStringSet(ddCENT_clean$nuc), 500)

n_worm <- length(worm_seq)
n_cent <- length(cent_seq)

#Step 3: Create pairwise combinations
#make indexes for the data, for iteration purposes
worm_idx <- rep(seq_along(worm_seq), each = n_cent)
cent_idx <- rep(seq_along(cent_seq), times = n_worm)

#Chunk size for each parallel task, to allow comparisons in parallel and to prevent RAM from being exceeded (and my computer exploding)
chunk_size <- 500  #starting chunk size, it worked so I didn't change it but it could be increased on computers with better RAM

#This line of code specifically is what does the splitting
chunks <- split(seq_along(worm_idx), ceiling(seq_along(worm_idx)/chunk_size))


#Step 4: Set up pairwise alignment to run in parallel
#Otherwise this will take forever to run

#This generalizes my code by checking for the number of cores available on the computer
num_cores <- parallel::detectCores() - 1 #it turns out I have 3 cores, who knew?
cl <- makeCluster(num_cores) #make a core cluster for parallel analysis
registerDoParallel(cl) #initiate parallel access


#Step 5: Parallel Pairwise Alignments
percent_id_list <- foreach(chunk = chunks, .combine = c, .packages = "Biostrings") %dopar% {
  Biostrings::pid(Biostrings::pairwiseAlignment(
    pattern = worm_seq[worm_idx[chunk]],
    subject = cent_seq[cent_idx[chunk]]
  ))
}

# Stop cluster
stopCluster(cl)

#Step 6: Convert the output into a matrix
similarity_matrix <- matrix(percent_id_list,
                            nrow = n_worm,
                            ncol = n_cent,
                            byrow = TRUE,
                            dimnames = list(paste0("W", seq_len(n_worm)),
                                            paste0("C", seq_len(n_cent))))

#Step 7: Prepare the data for the heatmap
heatmap_df <- melt(similarity_matrix)
colnames(heatmap_df) <- c("WORM_Index", "CENT_Index", "Percent_Identity")

#Step 8: Create the heatmap
heatmap_plot <- ggplot(heatmap_df, aes(x = WORM_Index, y = CENT_Index, fill = Percent_Identity)) +
  geom_tile() + #removed argument from geom_tile to get rid of borders around every alignment wasting valuable space
  scale_fill_gradient(low = "white", high = "purple", na.value = "grey90") +
  labs(
    title = "Pairwise % Nucleotide Identity Between Onychophora and Chilopoda",
    x = "Onychophora Sequences",
    y = "Chilopoda Sequences",
    fill = "% Identity"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    #axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text = element_blank(),      # removes all axis text (tick labels)
    axis.ticks = element_blank(),     # removes tick marks
    plot.title = element_text(face = "bold")
  )

# Save to JPEG
ggsave("similarity_heatmap.jpg", plot = heatmap_plot, width = 10, height = 8, dpi = 1200)
#I did this because the figure was too large to fit in the rstudio viewer. I also had to increase the dpi because the figure was blurry, but I don't think it made much of a difference



#1 Statistical analysis
# Simple t-test comparing average species count between the two groups
t_test_result <- t.test(species_count ~ group, data = combined_data)
print(t_test_result)
#This test adds a basic statistical comparison — showing whether Velvet Worms and Centipedes differ significantly in mean species richness per ecoregion

#provided a pictorial representation of the statistical report
# Add summary boxplot
ggplot(combined_data, aes(x = group, y = species_count, fill = group)) +
  geom_boxplot(alpha = 0.7) +
  scale_fill_manual(values = c("purple", "orange")) +
  labs(title = "Species Richness Comparison",
       subtitle = paste("p-value:", signif(t_test_result$p.value, 3)),
       x = "Group", y = "Species Count") +
  theme_minimal()

#2 Create a spatial heatmap to show areas where most records are concentrated.
#This helps identify whether both taxa overlap in ecological regions.
# Combine both datasets for comparison
all_data <- bind_rows(
  dfWORM_clean %>% mutate(group = "Onychophora"),
  ddCENT_clean %>% mutate(group = "Chilopoda")
)
coord_data <- all_data %>%
  separate(coord, into = c("Lat","Long"), sep = ",", convert = FALSE) %>%
  mutate(Lat  = as.numeric(str_replace_all(Lat, "[^0-9.-]", "")), Long = as.numeric(str_replace_all(Long, "[^0-9.-]", "")))
# Create a spatial heatmap showing density of records
ggplot(coord_data, aes(x = Long, y = Lat)) +
  geom_bin2d(bins = 50) +
  scale_fill_viridis_c(option = "D", direction = 1) +
  facet_wrap(~ group) +
  labs(
    title = "Species Occurrence Density by Geographic Location",
    subtitle = "Comparing Onychophora and Chilopoda Spatial Hotspots",
    x = "Longitude",
    y = "Latitude",
    fill = "Record Count"
  ) +
  theme_minimal(base_size = 12)

#3 To reduce repetition of codes, a function can take care of the two dataframes which makes the codes shorter
# This function filters data by ecoregion and species, 

# Create a reusable function to summarize species per group
summarize_species <- function(df, eco_shared, group_name) {
  df %>%
    filter(ecoregion %in% eco_shared) %>%
    group_by(ecoregion) %>%
    summarise(species_count = n_distinct(species)) %>%
    mutate(group = group_name)
}

# Apply the function to both datasets
worm_summary <- summarize_species(dfWORM, SharedEco, "Velvet Worm")
cent_summary <- summarize_species(ddCENT, SharedEco, "Centipede")
# Step 2: Combine into one dataframe
combined_data <- bind_rows(worm_summary, cent_summary)
combined_data <- combined_data %>% 
  mutate(ecoregion = gsub('_', ' ', ecoregion))
print(combined_data)

