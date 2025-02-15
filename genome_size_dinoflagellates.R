### How big are the genomes of dinoflagellates? ###
# V. POCHIC, 2025-02-15

# This little project was inspired by Tom White's post "Genome size" in his
# Datavision 2020 series (https://tom-e-white.com/datavision/). I downloaded 
# the NCBI data from his github 
# (https://github.com/tomwhite/datavision-code/tree/master/05-genome-size)
# and based my code on his. 
# Many thanks to him for making his excellent work public!

# Please remember to give credit if you're reusing this code (in parts or in 
# total)
# (CC-BY V. Pochic, 2025)

# Required packages ####
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(ggnewscale)

### Genomes from NCBI ####

# A function to extract information about organisms from the 'genomes' table
extract_group <- function(group_str) {
  groups <- str_split(group_str, ';')[[1]]
  if (groups[[1]] != "Eukaryota") {
    return(groups[[1]]);
  }
  if (groups[[2]] == "Fungi" | groups[[2]] == "Protists" | groups[[2]] == "Plants") {
    return(groups[[2]]);
  }
  return(groups[[3]]);
}

# Import the NCBI data
genomes <- read_csv("data/genomes.csv") %>%
  mutate(bases = `Size(Mb)` * 1000000) %>%
  rowwise() %>% # needed since following mutate uses a function
  mutate(group = extract_group(`Organism Groups`)) %>%
  filter(`Organism Groups` != "Eukaryota;Other;Other") %>%
  separate(`Organism Groups`, c("Organism1", "Organism2", "Organism3"), sep = ';') %>%
  filter(Organism3 != "Other Animals") %>%
  filter(bases > 0)

# Some exploration of the dataset
genomes %>% group_by(Organism1) %>% count()

genomes %>% group_by(group) %>% count()

counts <- genomes %>% group_by(Organism1, Organism2) %>% count()

### Dinoflagellates ####

# Import the data from Lin (2024)')'s table 1 : 'Published genomes of
# dinoflagellates'
dino_genomes_lin <- read.csv2('data/dino_genomes_lin_2024.csv', header = TRUE,
                              dec = '.', fileEncoding = 'ISO-8859-1')

# Modify the names to fit the other dataset
dino_genomes <- dino_genomes_lin %>%
  mutate(`#Organism Name` = Species) %>%
  # Transform genome size from Mbp to bp
  mutate(bases = Assembled.genome.size..Mbp.*10^6) %>%
  # Code variables that group organisms
  mutate(Organism1 = 'Eukaryota') %>%
  mutate(Organism2 = 'Protists') %>%
  # Code a variable indicating that it is a dinoflagellate
  mutate(dinoflagellate = 1) %>%
  # select only relevant variables
  select(`#Organism Name`, bases, Organism1, Organism2, dinoflagellate)

# In the other dataset
genomes2 <- genomes %>%
  # In this dataset, Chlorophytes are considered as Plants,
  # but we prefer to consider them as Protists
  mutate(Organism2 = ifelse(Organism3 == 'Green Algae',
                            'Protists',
                            Organism2)) %>%
  # And code the variable for dinoflagellates I only found 4 in the original
  # dataset (but maybe I missed some)
  mutate(dinoflagellate = ifelse(`#Organism Name` %in%
                                   c('Breviolum minutum', 'Cladocopium',
                                     'Symbiodinium',
                                     'Symbiodinium microadriaticum'),
                                 # These 3 are dinoflagellates
                                 1,
                                 # The other aren't
                                 0))%>%
  # Select only relevant variables
  select(`#Organism Name`, bases, Organism1, Organism2, dinoflagellate)

# Bind both datasets
genomes_combined <- bind_rows(genomes2, dino_genomes) %>%
  # We want to group animals together in 1 big 'Animals' category
  # The information we need for this is in the variable 'Organism2'
  mutate(Group = ifelse(Organism2 == 'Birds' |
                          Organism2 == 'Fishes' |
                          Organism2 == 'Mammals' |
                          Organism2 == 'Reptiles' |
                          Organism2 == 'Amphibians' |
                          Organism2 == 'Roundworms' |
                          Organism2 == 'Flatworms' |
                          Organism2 == 'Insects',
                        'Animals', 
                        # If it's not an eukaryota, it has to take the value 
                        # virus, bacteria, archaea, stored in 'Organism1'
                        ifelse(Organism1 != 'Eukaryota',
                               Organism1,
                               Organism2)
                        )
         ) %>%
  # Exclude an impossibly small animal genome indicated as Carcharodon carcharias 
  # (Great White Shark), with only 1005 bases
  filter(ifelse(Group == 'Animals' & bases < 2000, FALSE, TRUE)) %>%
  # Reorder 'Group' as we like
  mutate(Group = as_factor(Group)) %>%
  mutate(Group = fct_relevel(Group, c('Viruses', 'Archaea', 'Bacteria',
                                      'Fungi', 'Protists', 'Plants', 'Animals')))

# Let's check that
unique(genomes_combined$Group)
## Nice!

### Other organisms ####

# It could be nice to plot the genome size of other organisms, as a comparison.
# Here I chose some model organisms in biology and others that I find interesting,
# but there are many more!

# Get specific information on some organisms
models <- filter(genomes_combined, 
                 # Bacteria
                 # E. coli
                 (grepl('Escherichia coli', `#Organism Name`)
                  # Filter out E. coli phages
                  & Organism1 == 'Bacteria') |
                   # # Lactobacillus casei (found in cheese)
                   # (grepl('Lactobacillus casei', `#Organism Name`)
                   #  # Filter out phages
                   #  & Organism1 == 'Bacteria') |
                   # # Yersinia pestis (bubonic plague)
                   # (grepl('Yersinia pestis', `#Organism Name`)
                   #  # Filter out phages
                   #  & Organism1 == 'Bacteria')|
                   
                   ### Protists
                   # E. huxleyi (haptophyte)
                 (grepl('Emiliania huxleyi', `#Organism Name`)
                  # filter out phages
                  & Organism1 == 'Eukaryota') |
                   # # Plasmodium falciparum (malaria)
                   # grepl('Plasmodium falciparum', `#Organism Name`) |
                   # # Chlorella vulgaris (chlorophyte)
                   # grepl('Chlorella vulgaris', `#Organism Name`) |
                   
                   ### Plants
                   # Vitis vinifera (grape vine)
                   # (`#Organism Name` == 'Vitis vinifera') |
                   # # Triticum aestivum (wheat)
                   # grepl('Triticum aestivum', `#Organism Name`) |
                   # # Arabidopsis thaliana
                   grepl('Arabidopsis thaliana', `#Organism Name`) |
                   # Quercus suber (cork oak)
                   # grepl('Quercus suber', `#Organism Name`) |
                   
                   ### Animals
                   # Human
                   grepl('Homo sapiens', `#Organism Name`) |
                   # # Ambystoma mexicanum (axolotl)
                   # grepl('Ambystoma mexicanum', `#Organism Name`) |
                   # # Falco peregrinus (peregrine falcon)
                   # grepl('Falco peregrinus', `#Organism Name`) |
                   
                   ### Viruses
                   # Influenza
                   # grepl('Influenza B virus', `#Organism Name`) |
                   # # Rhinovirus (cold)
                   # grepl('Rhinovirus B14', `#Organism Name`) |
                   # # Nile crocodilepox virus (funny)
                   # grepl('Nile crocodilepox virus', `#Organism Name`) |
                   
                   ### Archaea
                   # Candidatus Heimdallarchaeota archaeon AB_125
                   # grepl('Candidatus Heimdallarchaeota archaeon AB_125', `#Organism Name`) |
                   
                   ### Fungi
                   # S. cerevisiae
                   (`#Organism Name` == 'Saccharomyces cerevisiae'
                    # filter out phages
                    & Organism2 == 'Fungi') |
                   # Tuber melanosporum (black truffle)
                   # (`#Organism Name` == 'Tuber melanosporum') |
                   
                   ### Get smallest and biggest dinoflagellate genomes
                   (bases == min(dino_genomes$bases) & dinoflagellate == 1) |
                   (bases == max(dino_genomes$bases) & dinoflagellate == 1)
                 
)




### Plotting ####

# Color palettes
# 7 groups, 7 colors
palette_7_groups <- c('#FF6448', '#791D40', '#977ED1', '#8D6456', '#11367D',
                      '#649003', '#FBA823')
# 2 group, 2 colors
palette_2_groups <- c('grey10', 'red3')
alpha_2_groups <- c(0.4,.8)

# The plot
ggplot() +
  
  ### Violin plots
  geom_violin(data = genomes_combined, aes(x = bases, y = Group, fill = Group,
                  color = Group), alpha = .5) +
  
  # custom color scale
  scale_color_discrete(type = palette_7_groups, guide = 'none') +
  scale_fill_discrete(type = palette_7_groups, guide = 'none') +
  
  ### Dinoflagellates
  # indicate specifically dinoflagellates with red dots
  geom_point(data = subset(genomes_combined, dinoflagellate == 1),
             aes(x = bases, y = Group), color = 'red3', size = 1.5,
             shape = 21, fill = 'red3') +
  
  ### Other organisms
  # New color scale
  new_scale_color() +
  scale_color_discrete(type = palette_2_groups, guide = 'none') +
  
  # Dots that indicate genome size (for our organisms in 'models')
  geom_point(data = subset(models, dinoflagellate == 0),
             aes(x = bases, y = Group), color = 'grey10', size = 1,
             shape = 21, fill = 'grey10') +
  # Labels associated with each data point
  geom_label_repel(data = models,
            aes(x = bases, y = Group, label = `#Organism Name`,
                color = as_factor(dinoflagellate),
                alpha = as_factor(dinoflagellate)
                ),
            size = 2.5, min.segment.length = 0,
            force_pull = 0.1) +
  
  # alpha scale
  scale_alpha_manual(values = alpha_2_groups, guide = 'none') +
  
  # log scale on y-axis
  scale_x_log10(
    breaks = c(10^2, 10^3, 10^4, 10^5, 10^6, 10^7, 10^8, 10^9, 10^10),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  # Add a manmade legend :)
  geom_point(aes(x = 10^9, y = 'Viruses'), 
             shape = 21, color = 'red3', fill = 'red3',
             size = 2,
             stroke = .8) +
  geom_text(aes(x = 10^9.8, y = 'Viruses', label = 'Dinoflagellates'), size = 2.8,
            color = 'red3') +
  # Labels
  labs(
    title = "Genome size (in base pairs)", x = NULL, y = NULL,
    caption = "Data source: 48790 genomes on NCBI (2019) and 30 dinoflagellate genomes in Lin (2024)"
  ) +
  # Theme
  theme_classic() + 
  theme(axis.title.y = element_blank(), 
        panel.grid.major.x = element_line(linewidth = .1, color = 'grey60'))

# NB : Because the function geom_label_repel() contains a random element, the
# positions of the labels vary slightly each time you run the plot function

# Save as an image
# ggsave("genome-size_dinos.tiff", width = 164, height = 148,
#        units = 'mm', dpi = 600, compression = 'lzw')
