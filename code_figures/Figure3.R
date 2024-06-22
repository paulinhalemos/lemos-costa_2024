# Figure 3 - pvalues

dt <- read_csv("table_results.csv") %>% filter(model == 1) %>% 
  mutate(pv = paste(row_number(), model, collapse = "-")) %>% mutate(ds = substr(label, 1, 3)) %>% 
  select(pv, ds, pval_rnd, pval_scramble) %>% pivot_longer(names_to = "type", values_to = "pvalue", cols = -c(pv, ds))

dt <- dt %>% mutate(ds = ifelse(ds == "Cad", "Cadotte", ifelse(ds == "Bio", "Biodiversity II", "Wageningen"))) %>% rename(`data set` = ds)

ggplot(dt) + aes(x = pvalue, fill = `data set`) +  geom_histogram(breaks = seq(0, 0.5, by = 0.05), color = "white") + theme_bw() + theme(legend.position = "bottom") + 
  scale_fill_manual(values = c("#8ecae6", "#219ebc", "#023047"))
