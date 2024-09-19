###------------------------------------------------------###
###----------- Exploratory Data Analysis  ---------------###
###---------- (Supplementary Section D.1)  --------------###
###------------------------------------------------------###

data = read.csv("cholera_data.csv")

###--- Bari size distribution ---###
Ni.bari = data %>% 
  group_by(id) %>%
  summarise(Ni = n())

nrow(Ni.bari)  ## Number of baris

table(Ni.bari$Ni)

ggplot(Ni.bari, aes(x = Ni)) +
  geom_histogram(breaks = seq(1, 250, by = 3), 
                 fill = "lightblue", 
                 color = "black",
                 size = 0.1) +
  labs(#title = "Histogram of Bari (Household) size (Total: 5625 Baris)",
    x = "Bari size", 
    y = "Number of Baris") +
  annotate("text", 
           x = 239, y = 100, 
           label = "Bari size 239 \n \u2193", 
           hjust = 0.5) + 
  theme_bw()

# Save the plot to a PDF file
ggsave("FigS5.Bari_size_histogram.jpeg", width = 8, height = 2, units = "in")


###--- Observed time distribution ---###
data = data %>% mutate(T = ifelse(D == 1, Y, NA), C = ifelse(D == 0, Y, NA))

p.Y.dist = ggplot(data, aes(x = Y)) +
  geom_histogram(bins = 100,
                 fill = "lightgreen", 
                 color = "black") +
  labs(#title = "Histogram of Observed time", 
    x = "Days", 
    y = "Frequency") +
  theme_minimal()

###--- Event time distribution ---###
p.T.dist = ggplot(data, aes(x = T)) +
  geom_histogram(bins = 100,
                 fill = "lightblue", 
                 color = "black",
                 size = 0.1) +
  xlim(c(0, 470)) +
  labs(#title = "Histogram of Observed time", 
    x = "Event Time (Days)", 
    y = "Frequency") +
  theme_minimal()

###--- Censoring time distribution ---###
p.C.dist = ggplot(data, aes(x = C)) +
  geom_histogram(bins = 100,
                 fill = "lightpink", 
                 color = "black",
                 size = 0.1) +
  xlim(c(0, 470)) +
  labs(#title = "Histogram of Observed time", 
    x = "Censoring Time (Days)", 
    y = "Frequency") +
  theme_minimal()

plot_grid(p.T.dist, p.C.dist, ncol = 1, align = "v")
ggsave("FigS6.Event_time_censoring_time_distribution.pdf", width = 8, height = 6)


###--- Vaccination coverage distribution ---###
ggplot(data %>%
         group_by(id) %>%
         summarise(A.bar = mean(A)), 
       aes(x = A.bar)) +
  geom_histogram(bins = 101,
                 fill = "lightgreen", 
                 color = "black",
                 size = 0.1) +
  labs(x = "Vaccination coverage", 
       y = "Number of Baris") +
  theme_bw()

# Save the plot to a PDF file
ggsave("FigS7.Vaccine_coverage_histogram.pdf", width = 8, height = 2)