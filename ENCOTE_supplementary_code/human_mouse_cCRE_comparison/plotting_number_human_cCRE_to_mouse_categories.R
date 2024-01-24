#R 3.5.1
library(ggplot2)
library(scales)
source('reference_files/ggplotting_source', encoding = 'UTF-8')

#Make bar charts for total human-mouse cCRE categories
Categories <- c("Shared\ncCRE", "Different\ncCRE", "Orthologous TE",
                "Syntenic\nHuman cCRE", "Non-syntenic\nHuman cCRE", "Human\nTE-cCRE")
Numbers <- c(96792, 34237, 18010, 394610, 167134, 215752)
ccre_data <- data.frame(Categories, Numbers)
ccre_data <- within(ccre_data, Categories <- factor(Categories,
                                                    levels=Categories))
#Plot
ggplot(ccre_data, aes(Categories)) +
  geom_bar(aes(weight = Numbers)) +
  scale_y_continuous(labels=scales::label_number(scale=0.001, suffix="K")) +
  labs(y='Number', x='')

#Combine Shared/Different cCRE categories
Categories <- c("Syntenic\ncCRE", "Orthologous TE",
                "Syntenic\nHuman cCRE", "Non-syntenic\nHuman cCRE", "Human\nTE-cCRE")
Numbers <- c(131029, 18010, 394610, 167134, 215752)
ccre_data <- data.frame(Categories, Numbers)
ccre_data <- within(ccre_data, Categories <- factor(Categories,
                                                    levels=Categories))
#Plot (saved as 4x10 for small, 6x10 for normal)
ggplot(ccre_data, aes(Categories)) +
  geom_bar(aes(weight = Numbers)) +
  scale_y_continuous(labels=scales::label_number(scale=0.001, suffix="K")) +
  geom_text(aes(label=Numbers, y=Numbers), vjust=-0.25) +
  labs(y='Number', x='')

