library(tidyverse)
library(ggdist)
library(ggsvg)

dbRES <- read_csv("~/Desktop/dbRES.csv")

plot_data <- dbRES %>% 
  group_by(phylum) %>%
  mutate(median_res = median(res, na.rm=T),
         n = n()) %>%
  ungroup() %>%
  arrange(desc(kingdom),median_res,phylum) %>%
  mutate(phylum = factor(phylum, levels=unique(.$phylum))) %>% 
  mutate(image = case_when(phylum == "Acanthocephala" ~ "https://images.phylopic.org/images/4f117460-d00a-49e4-865f-52bcc7d5887b/vector.svg",
                         phylum == "Annelida" ~         "https://images.phylopic.org/images/079b4cee-ba72-4105-b59c-59c5c9591549/vector.svg",
                         phylum == "Arthropoda" ~       "https://images.phylopic.org/images/f9b57bb5-7675-4c53-936f-924e8f5a75f3/vector.svg",
                         phylum == "Chordata" ~         "https://images.phylopic.org/images/2ca942d6-57cd-4281-89c2-42350491984f/vector.svg",
                         phylum == "Mollusca" ~         "https://images.phylopic.org/images/443c7494-aa78-4a21-b2b4-cfa3639e1346/vector.svg",
                         phylum == "Platyhelminthes" ~  "https://images.phylopic.org/images/b81d728e-0e8c-4faf-8f40-edf999143f10/vector.svg",
                         phylum == "Porifera" ~         "https://images.phylopic.org/images/3449d9ef-2900-4309-bf22-5262c909344b/vector.svg",
                         phylum == "Rotifera"~          "https://images.phylopic.org/images/3042a73d-353d-4191-811f-9b12f57c958c/vector.svg",
                         phylum == "Tardigrada" ~       "https://images.phylopic.org/images/0957b8c0-a052-4f7d-a1ac-95f2179c6582/vector.svg",
                         phylum == "Rhodophyta" ~       "https://images.phylopic.org/images/e9df48fe-68ea-419e-b9df-441e0b208335/vector.svg",
                         phylum ==  "Tracheophyta" ~    "https://images.phylopic.org/images/585f5a9d-b96d-4fdc-aec7-01b9540e057f/vector.svg",
                         phylum == "Nematoda" ~         "https://images.phylopic.org/images/3c60fbfb-5722-4248-94f0-23f841030294/vector.svg" ,
                         phylum == "Bryophyta" ~        "https://images.phylopic.org/images/c9433947-a86d-453e-817a-7aba27fb453f/vector.svg" , 
                         phylum == "Brachiopoda" ~      "https://images.phylopic.org/images/c0daf6f8-876d-4b61-98ad-fc8cde084841/vector.svg",
                         phylum == "Bryozoa" ~          "https://images.phylopic.org/images/191ad6ce-78f2-444d-b4ad-9c4162b1803f/vector.svg" ,
                         phylum == "Chaetognatha" ~     "https://images.phylopic.org/images/345d2e9d-c443-4579-8f5b-c9b15ab55c84/vector.svg",
                         phylum == "Cnidaria" ~         "https://images.phylopic.org/images/26c0169f-b4a2-4871-8b30-e00db6e5958d/vector.svg" ,
                         phylum == "Ctenophora" ~       "https://images.phylopic.org/images/2fa866ea-fa23-4b22-9382-66139a9c2cf1/vector.svg" ,
                         phylum == "Cycliophora" ~      "https://images.phylopic.org/images/d08250be-c00c-48da-a76f-77eee31b09dc/vector.svg",
                         phylum == "Echinodermata" ~    "https://images.phylopic.org/images/77fa45d8-a9e5-4067-b2df-ffe4c7f0cab2/vector.svg" ,
                         phylum == "Hemichordata" ~     "https://images.phylopic.org/images/369c9a71-fea9-4b5c-8718-07735eaba5fc/vector.svg",
                         phylum == "Kinorhyncha" ~      "https://images.phylopic.org/images/bc37a697-84bf-4e22-b919-de09657b4e93/vector.svg",
                         phylum == "Loricifera" ~       "https://images.phylopic.org/images/03219505-5777-42e6-9660-9ddb80746144/vector.svg",
                         phylum == "Nematomorpha" ~     "https://images.phylopic.org/images/23c6cb9b-855e-4dd2-8f24-d8882a24db30/vector.svg",
                         phylum == "Nemertea" ~         "https://images.phylopic.org/images/1ae2caed-ef2a-4151-a023-17b8c601d671/vector.svg",
                         phylum == "Priapulida" ~       "https://images.phylopic.org/images/4c501503-6407-4534-8cee-0d33aa2e6fbd/vector.svg" ,
                         phylum == "Basidiomycota" ~    "https://images.phylopic.org/images/0ec7e6e6-e115-450a-b346-c40685ed3b1a/vector.svg",
                         phylum == "Anthocerotophyta" ~ "https://images.phylopic.org/images/23c6cb9b-855e-4dd2-8f24-d8882a24db30/vector.svg",
                         phylum == "Marchantiophyta" ~  "https://images.phylopic.org/images/4054a966-b5e1-425d-924e-352593ff9a1d/vector.svg")
  ) %>% 
  mutate(pic.x.pos = 9.5)
   

plot_data %>%
  ggplot(aes(
     x = res,
     y = phylum,
     fill = kingdom,
     color = kingdom
  )) +
  geom_boxplot(
    data = ~ .x %>%  filter(n == 1),
    aes(x = median_res, fill = kingdom),
    colour = "black",
  ) +
  stat_slab(
    data = ~ .x %>% filter(n > 1),
    color = "gray15",
    size = .3,
    expand = FALSE,
    trim = TRUE,
    height = 3
  ) +
  geom_point_svg(data = ~ .x %>% 
                    distinct(phylum, .keep_all=TRUE) %>% 
                    rowwise() %>% 
                    mutate(svg = paste(readLines(as.character(image)), collapse="\n")) %>% 
                    ungroup(),
                    aes(x=jitter(pic.x.pos, factor=.00001),y = phylum, svg = I(svg))) +
  labs(x = "Residuals", y = NULL) +
  xlim(-10, 10)+
  annotate(
    "segment",
    x = 0.5,
    xend = 3.5,
    y = 30,
    yend = 30,
    color = "grey30",
    arrow = arrow(
      ends = "last",
      angle = 15,
      length = unit(.2, "cm")
    )
  ) +
  annotate(
    "text",
    x = 0.5,
    y = 30.5,
    hjust = 0,
    vjust = 0.5,
    size = 5,
    color = "grey30",
    label = "Popular interest"
  ) +
  
  annotate(
    "segment",
    x = -0.5,
    xend = -3.5,
    y = 30,
    yend = 30,
    color = "grey30",
    arrow = arrow(
      ends = "last",
      angle = 15,
      length = unit(.2, "cm")
    )
  ) +
  annotate(
    "text",
    x = -0.5,
    y = 30.5,
    hjust = 1,
    vjust = 0.5,
    size = 5,
    color = "grey30",
    label = "Scientific interest"
  ) +
  
  # geom_vline(lty = 3, size = 0.5, col = "black", xintercept = 0) +
  # scale_color_manual(values = custom_color)+
  # scale_fill_manual(values = custom_color)+
  scale_fill_viridis_d() +
  scale_y_discrete(drop = FALSE, expand = c(0.05, .05)) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(size = 12)
  )
