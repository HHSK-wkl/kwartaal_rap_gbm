# Libraries
library(tidyverse)
library(knitr)
library(DT)
library(leaflet)
library(HHSKwkl)
library(readxl)
library(lubridate)
library(gt)
library(glue)
library(ggbeeswarm)
library(scales)
# library(tsibble)

# Other options

options(OutDec = ",")

# ggplot2

theme_set(hhskthema())

# DT

dt_labels_nederlands()

my_datatable <- function(df, ...) {
  DT::datatable(data = df, extensions = 'Buttons',
            options = list(dom = 'lfirtpB', buttons = c('csv', 'excel', 'pdf')), ...)
}

## ---- load-data ----

q_eind <- ymd(20250101)
q_begin <- q_eind - period(3, "month")
q_tekst <- as.character(tsibble::yearquarter(q_begin))
# q_tekst <- "2024 Q2+Q3"

# copy_data(c("fys_chem.rds", "meetpunten.csv", "parameters.csv", "normen.txt", "gbm_toxiciteit.xlsx"))
# download_data(c("fys_chem.rds", "meetpunten.csv", "parameters.csv", "normen.txt"))
# HHSKwkl::download_data(c("fys_chem.rds", "meetpunten.rds", "parameters.rds", "normen.rds", "gbm_toelating_werking.rds", "gbm_middelen.rds"))
fys_chem <- readRDS("data/fys_chem.rds")
meetpunten <-readRDS("data/meetpunten.rds")
parameters <- readRDS("data/parameters.rds")
normen <- readRDS("data/normen.rds") %>%
  group_by(parnr) %>%
  mutate(norm = min(norm_JGM, norm_MAX, norm_P90, na.rm = TRUE),
         normtype = case_when(norm == norm_JGM ~ "JG-MKN", norm == norm_MAX ~ "MAC-MKN", norm == norm_P90 ~ "MTR"))
toxiciteit <- read_excel("data/gbm_toxiciteit.xlsx", sheet = "SSDinfo")
toelatingen <- readRDS("data/gbm_toelating_werking.rds")
middelen <- readRDS("data/gbm_middelen.rds")
  

# ws_grens <- sf::st_read("data/ws_grens.gpkg", crs = 28992)

gbm_meetpunten <- 
  meetpunten %>% 
  filter(landgebruik == "glastuinbouw" | str_detect(mp, "^GGA")) %>% 
  filter(!mp %in% c("S_0683"))
  
# gbm_meetpunten <- meetpunten %>% filter(c12 == 2 | str_detect(mp, "^GGA"))

f_aquopar <- maak_opzoeker(parameters, parnr, aquo_parcode)

opzoek_pars <- parameters %>%
  select(parnr, parnaamlang) %>%
  mutate(parnaamlang = str_to_sentence(parnaamlang))%>%
  tibble::deframe()


mspaf <- function(paf_vector){
  1 - prod(1 - paf_vector, na.rm = TRUE)
}

data_gbm <-
  fys_chem %>%
  semi_join(gbm_meetpunten, by = "mp") %>%
  filter(parnr >= 1000, parnr <= 1999, 
         year(datum) >= year(Sys.Date()) - 3,
         datum < q_eind) %>% 
  mutate(kwartaal = tsibble::yearquarter(datum),
         kwartaal_c = as.character(kwartaal)) %>%
  # filter(kwartaal <= tsibble::yearquarter(today())) %>%
  mutate(laatst = (datum >= q_begin & datum < q_eind)) %>%
  mutate(datum = case_when(
    datum == ymd("20221024") ~ ymd("20221101"),
    TRUE ~ datum
  ))


data_n <- fys_chem %>%
  add_jaar_maand() %>%
  filter(parnr == 7, jaar > 2010) %>%
  mutate(kwartaal = quarter(datum)) %>%
  mutate(jaarmaand = tsibble::yearmonth(datum)) %>%
  left_join(meetpunten, by = "mp") %>%
  mutate(glas = landgebruik == "glastuinbouw") %>%
  filter(!is.na(glas)) %>%
  filter(datum < q_eind) %>%
  mutate(tekst = ifelse(glas, "Glastuinbouw", "Overige gebieden"))

data_p <- fys_chem %>%
  add_jaar_maand() %>%
  filter(parnr == 3, jaar > 2010) %>%
  mutate(kwartaal = quarter(datum)) %>%
  mutate(jaarmaand = tsibble::yearmonth(datum)) %>%
  left_join(meetpunten, by = "mp") %>%
  mutate(glas = landgebruik == "glastuinbouw") %>%
  filter(!is.na(glas)) %>%
  filter(datum < q_eind) %>%
  mutate(tekst = ifelse(glas, "Glastuinbouw", "Overige gebieden"))

## ---- maak-figuren ----

plot_mspaf_acuut <-
  data_gbm %>%
  mutate(paf_acuut = paf_gbm(f_aquopar(parnr),
                             concentratie = waarde,
                             detectiegrens = detectiegrens,
                             ssd_data = toxiciteit,
                             type_paf = "acuut"),
         paf_chronisch = paf_gbm(f_aquopar(parnr),
                                 concentratie = waarde,
                                 detectiegrens = detectiegrens,
                                 ssd_data = toxiciteit,
                                 type_paf = "chronisch")) %>%
  group_by(mp, datum, kwartaal, kwartaal_c, laatst) %>%
  summarise(mspaf_acuut = mspaf(paf_acuut),
            mspaf_chronisch = mspaf(paf_chronisch)) %>%
  mutate(mspaf_acuut_i = 1 / mspaf_acuut) %>%
  filter(mspaf_acuut_i < 10000) %>%

  ggplot(aes(kwartaal_c, mspaf_acuut_i)) +
  geom_hline(yintercept = c(10, 200), colour = oranje, linetype = c(1,2)) +
  geom_quasirandom(colour = "grey50", size = 1.8, width = 0.25) +
  scale_y_continuous(transform = c("log10", "reverse"), breaks = breaks_log(n = 8), 
                     labels = scales::label_number(prefix = "1 op de ", big.mark = ".", accuracy = 1),
                     sec.axis = sec_axis(transform = ~ 1 / ., name = "msPAF acuut", breaks = breaks_log(n = 8),
                                         labels = function(x) scales::percent(x, decimal.mark = ",", accuracy = 0.001))) +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  labs(title = "Potentieel aangetaste fractie soorten - acuut (msPAF)",
       subtitle = "per monster per kwartaal",
       y = "Aandeel aangetaste soorten",
       x = "",
       caption = "Monsters waar minder dan 1 op de 10.000 soorten worden aangetast zijn niet getoond.") +
  theme(plot.subtitle = element_text(face = "plain"),
        panel.grid.major.x = element_blank(),
        axis.ticks.x = element_blank()) +
  guides(fill = "none") +
  NULL

plot_mspaf_chronisch <-
  data_gbm %>%
  mutate(paf_acuut = paf_gbm(f_aquopar(parnr),
                             concentratie = waarde,
                             detectiegrens = detectiegrens,
                             ssd_data = toxiciteit,
                             type_paf = "acuut"),
         paf_chronisch = paf_gbm(f_aquopar(parnr),
                                 concentratie = waarde,
                                 detectiegrens = detectiegrens,
                                 ssd_data = toxiciteit,
                                 type_paf = "chronisch")) %>%
  group_by(mp, datum, kwartaal, kwartaal_c, laatst) %>%
  summarise(mspaf_acuut = mspaf(paf_acuut),
            mspaf_chronisch = mspaf(paf_chronisch)) %>%
  mutate(mspaf_chronisch_i = 1 / mspaf_chronisch) %>%
  filter(mspaf_chronisch_i < 10000) %>%
  ggplot(aes(kwartaal_c, mspaf_chronisch_i)) +
  geom_hline(yintercept = c(20, 200), colour = oranje, linetype = c(1,2)) +
  geom_quasirandom(colour = "grey50", size = 1.8, width = 0.25) +
  scale_y_continuous(transform = c("log10", "reverse"), breaks = breaks_log(n = 8), 
                     labels = scales::label_number(prefix = "1 op de ", big.mark = ".", accuracy = 1),
                     sec.axis = sec_axis(transform = ~ 1 / ., name = "msPAF chronisch", breaks = breaks_log(n = 8),
                                         labels = function(x) scales::percent(x, decimal.mark = ",", accuracy = 0.001))) +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  labs(title = "Potentieel aangetaste fractie soorten - chronisch (msPAF)",
       subtitle = "per monster per kwartaal",
       y = "Aandeel aangetaste soorten",
       x = "",
       caption = "Monsters waar minder dan 1 op de 10.000 soorten worden aangetast zijn niet getoond.") +
  theme(plot.subtitle = element_text(face = "plain"),
        panel.grid.major.x = element_blank(),
        axis.ticks.x = element_blank()) +
  guides(fill = "none") +
  NULL


plot_perc_norm <-
  data_gbm %>%
  left_join(normen, by = "parnr") %>%
  mutate(overschrijding = waarde > norm & is.na(detectiegrens)) %>%
  group_by(kwartaal_c, laatst) %>%
  summarise(n = n(),
            n_ov = sum(overschrijding, na.rm = TRUE),
            frac = n_ov / n()) %>%
  ggplot(aes(kwartaal_c, frac, fill = laatst)) +
  geom_col() +
  geom_text(aes(label = glue("{n_ov}
                             van
                             {n}")), y = 0.0005, colour = "white", fontface = "bold", vjust = 0) +
  scale_fill_manual(values = c("grey60", hhskblauw)) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(c(0, 0.1)),
                     labels = function(x) scales::percent(x, accuracy = 0.10, decimal.mark = ",")) +
  # scale_x_date(date_breaks = "year", date_labels = "%Y") +
  labs(title = "Gewasbeschermingsmiddelenmetingen boven norm",
       subtitle = "Percentage van alle metingen",
       y = "",
       x = "") +
  theme(plot.subtitle = element_text(face = "plain"),
        panel.grid.major.x = element_blank(),
        axis.ticks.x = element_blank()) +
  guides(fill = FALSE)

plot_aantal_norm <-
  data_gbm %>%
  left_join(normen, by = "parnr") %>%
  filter(waarde > norm & is.na(detectiegrens)) %>%
  group_by(kwartaal_c, laatst) %>%
  summarise(overschrijdende_stoffen = n_distinct(parnr)) %>%

  ggplot(aes(kwartaal_c, overschrijdende_stoffen, fill = laatst)) +
  geom_col() +
  geom_text(aes(label = overschrijdende_stoffen), y = 0.5, colour = "white", fontface = "bold", vjust = 0) +
  scale_fill_manual(values = c("grey60", hhskblauw)) +
  scale_y_continuous(limits = c(0,NA), expand = expansion(c(0, 0.1))) +
  # scale_x_date(date_breaks = "year", date_labels = "%Y") +
  labs(title = "Aantal gewasbeschermingsmiddelen boven norm",
       subtitle = "per kwartaal",
       y = "Aantal stoffen",
       x = "") +
  theme(plot.subtitle = element_text(face = "plain"),
        panel.grid.major.x = element_blank(),
        axis.ticks.x = element_blank()) +
  guides(fill = FALSE)



plot_sno <-
  data_gbm %>%
  left_join(normen, by = "parnr") %>%
  filter(is.na(detectiegrens)) %>%
  mutate(overschrijdingsfactor = waarde / norm) %>%
  group_by(kwartaal_c, laatst, mp, datum) %>%
  filter(overschrijdingsfactor > 1) %>%
  summarise(som_factor = sum(overschrijdingsfactor, na.rm = TRUE)) %>%
  group_by(kwartaal_c, laatst) %>%
  summarise(gem_factor = mean(som_factor, na.rm = TRUE)) %>%

  ggplot(aes(kwartaal_c, gem_factor, fill = laatst)) +
  geom_col() +
  geom_text(aes(label = round(gem_factor)), y = 0.1, colour = "white", fontface = "bold", vjust = 0) +
  scale_fill_manual(values = c("grey60", hhskblauw)) +
  # scale_y_log10(expand = expansion(mult = c(0, 0.1))) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(c(0, 0.1))) +
  labs(title = "Gemiddelde totale overschrijdingsfactor",
       subtitle = "per kwartaal - op logaritmische schaal",
       caption = "De totale overschrijdingsfactor is de som van alle overschrijdingsfactoren in een monster",
       y = "Overschrijdingsfactor",
       x = "") +
  theme(plot.subtitle = element_text(face = "plain"),
        panel.grid.major.x = element_blank(),
        axis.ticks.x = element_blank()) +
  guides(fill = FALSE) +
  # geom_vline(xintercept = ymd(20200101, 20190101, 20180101), colour = "grey70", linetype = "dashed") +
  NULL

plot_aantal <-
  data_gbm %>%
  group_by(kwartaal_c, laatst) %>%
  filter(is.na(detectiegrens)) %>%
  summarise(n_stoffen = n_distinct(parnr)) %>%

  ggplot(aes(kwartaal_c, n_stoffen, fill = laatst)) +
  geom_col() +
  geom_text(aes(label = n_stoffen), y = 1, colour = "white", fontface = "bold", vjust = 0) +
  scale_fill_manual(values = c("grey60", hhskblauw)) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(c(0, 0.1))) +
  # scale_x_date(date_breaks = "year", date_labels = "%Y") +
  labs(title = "Aantal aangetroffen gewasbeschermingsmiddelen",
       subtitle = "verschillende stoffen per kwartaal",
       y = "Aantal verschillende stoffen",
       x = "") +
  theme(plot.subtitle = element_text(face = "plain"),
        panel.grid.major.x = element_blank(),
        axis.ticks.x = element_blank()) +
  guides(fill = FALSE)


plot_monster <-
  data_gbm %>%
  filter(is.na(detectiegrens)) %>%
  group_by(kwartaal_c, laatst, mp, datum) %>%
  summarise(n_monster = n_distinct(parnr)) %>%
  group_by(kwartaal_c, laatst) %>%
  summarise(gem_per_monster = mean(n_monster)) %>%

  ggplot(aes(kwartaal_c, gem_per_monster, fill = laatst)) +
  geom_col() +
  geom_text(aes(label = round(gem_per_monster, 1)), y = 0.5, colour = "white", fontface = "bold", vjust = 0) +
  scale_fill_manual(values = c("grey60", hhskblauw)) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(c(0, 0.1))) +
  # scale_x_date(date_breaks = "year", date_labels = "%Y") +
  labs(title = "Aantal gewasbeschermingsmiddelen per monster",
       subtitle = "gemiddeld per kwartaal",
       y = "Aantal stoffen",
       x = "") +
  theme(plot.subtitle = element_text(face = "plain"),
        panel.grid.major.x = element_blank(),
        axis.ticks.x = element_blank()) +
  guides(fill = FALSE)

# tabel_overschr <-
#   data_gbm %>%
#   left_join(normen, by = "parnr") %>%
#   filter(laatst,
#          is.na(detectiegrens),
#          waarde > norm) %>%
#   mutate(overschrijdingsfactor = signif(waarde / norm, digits = 2),
#          parnaam = opzoek_pars[as.character(parnr)]) %>%
#   arrange(desc(overschrijdingsfactor)) %>%
#   select(Meetpunt = mp,
#          Datum = datum,
#          Normtype = normtype,
#          Stofnaam = parnaam,
#          Overschrijdingsfactor = overschrijdingsfactor)
#
# tabel_overschr %>%
#   group_by(Stofnaam)  %>%
#   summarise(`Aantal keer boven norm` = n(),
#             `Gemiddelde overschrijdingsfactor` = signif(mean(Overschrijdingsfactor), digits = 2)) %>%
#   arrange(desc(`Aantal keer boven norm`))

tabel_basis <-
  data_gbm %>%
  left_join(normen, by = "parnr") %>%
  filter(laatst, is.na(detectiegrens)) %>%
  mutate(parnaam = opzoek_pars[as.character(parnr)],
         overschrijding = waarde > norm,
         overschrijdingsfactor = waarde / norm,
         overschrijdingsfactor = ifelse(is.na(overschrijdingsfactor) | overschrijdingsfactor <= 1, 0, overschrijdingsfactor)) %>%
  mutate(paf_acuut = paf_gbm(f_aquopar(parnr),
                             concentratie = waarde,
                             detectiegrens = detectiegrens,
                             ssd_data = toxiciteit,
                             type_paf = "acuut"),
         paf_chronisch = paf_gbm(f_aquopar(parnr),
                                 concentratie = waarde,
                                 detectiegrens = detectiegrens,
                                 ssd_data = toxiciteit,
                                 type_paf = "chronisch")) %>%
  left_join(toelatingen) %>% 
  mutate(toegelaten = if_else(toegelaten, "Ja", "Nee", missing = "Nee")) %>% 

  group_by(parnr, parnaam, toegelaten, werking) %>%
  summarise(n_aanwezig = n(),
            n_norm = sum(overschrijding, na.rm = TRUE),
            n_tox = sum(paf_acuut > 0.005 | paf_chronisch > 0.005, na.rm = TRUE),
            gem_factor = sum(overschrijdingsfactor * overschrijding) / n_norm,
            max_factor = max(overschrijdingsfactor * overschrijding),
            max_paf_acuut = max(paf_acuut),
            max_paf_chronisch = max(paf_chronisch)) %>%
  ungroup() %>% 
  mutate(gem_factor = ifelse(is.na(gem_factor), 0, round(gem_factor, digits = 1)),
         max_factor = ifelse(is.na(max_factor), 0, round(max_factor, digits = 1))) %>%
  arrange(desc(n_norm), desc(gem_factor)) %>%
  filter(n_norm > 0 | n_tox > 0)

tabel <- tabel_basis %>% 
  select(Stofnaam = parnaam,
         `Aantal metingen boven norm` = n_norm,
         `Aantal metingen toxisch (> 0,5%)` = n_tox,
         `Gemiddelde overschrijding` = gem_factor,
         `Hoogste overschrijding` = max_factor,
         `Hoogste PAF (acuut)` = max_paf_acuut,
         `Hoogste PAF (chronisch)` = max_paf_chronisch,
         Toegelaten = toegelaten,
         Werking = werking) # , `Aantal keer aangetroffen` = n_aanwezig


tabel_meetwaarden <-
  data_gbm %>%
  left_join(normen, by = "parnr") %>%
  filter(laatst, is.na(detectiegrens)) %>%
  mutate(parnaam = opzoek_pars[as.character(parnr)],
         overschrijding = waarde > norm,
         overschrijdingsfactor = waarde / norm
         ) %>%
  mutate(paf_acuut = paf_gbm(f_aquopar(parnr),
                             concentratie = waarde,
                             detectiegrens = detectiegrens,
                             ssd_data = toxiciteit,
                             type_paf = "acuut"),
         paf_chronisch = paf_gbm(f_aquopar(parnr),
                                 concentratie = waarde,
                                 detectiegrens = detectiegrens,
                                 ssd_data = toxiciteit,
                                 type_paf = "chronisch")) %>%
  filter(overschrijding | paf_acuut > 0.005 | paf_chronisch > 0.005) %>%
  arrange(desc(paf_acuut)) %>%
  # mutate(paf_acuut = scales::percent(paf_acuut, decimal.mark = ",", accuracy = 0.1),
         # paf_chronisch = scales::percent(paf_chronisch, decimal.mark = ",", accuracy = 0.1)) %>%
  select(Stofnaam = parnaam,
         Meetlocatie = mp,
         Datum = datum,
         `Meetwaarde (ug/l)` = waarde,
         Overschrijdingsfactor = overschrijdingsfactor,
         `PAF acuut` = paf_acuut,
         `PAF chronisch` = paf_chronisch) %>%
  DT::datatable(extensions = 'Buttons', rownames = FALSE, filter = "top",
                options = list(dom = 'lfirtpB', buttons = c('csv', 'excel', 'pdf'), pageLength = 50)) %>%
  DT::formatRound(c(4), dec.mark = ",", digits = 3) %>%
  DT::formatRound(c(5), dec.mark = ",", digits = 1, mark = "" ) %>%
  DT::formatPercentage(c(6, 7), dec.mark = ",", digits = 1)

tabel_middelen <- 
  middelen %>% 
  semi_join(tabel_basis, by = "parnr") %>% 
  mutate(Toegelaten = if_else(expiratiedatum > q_eind, "Ja", "Nee", missing = "Nee")) %>% 
  arrange(parnaamlang, desc(expiratiedatum)) %>% 
  select(Stofnaam = parnaamlang,
         Middelnaam = middelnaam,
         `Werkzame stoffen` = werkzame_stoffen_middel,
         Expiratiedatum = expiratiedatum,
         Toegelaten,
         Werking = aard_werking_middel
          ) %>% 
  DT::datatable(extensions = 'Buttons', rownames = FALSE, filter = "top",
                options = list(dom = 'lirtpB', buttons = c('csv', 'excel', 'pdf'), pageLength = 50))

# data_n %>%
#   group_by(kwartaal_c, laatst) %>%
#   summarise(gem_n = mean(waarde)) %>%
#   ggplot(aes(kwartaal_c, gem_n, fill = laatst)) +
#   geom_col() +
#
#   scale_fill_manual(values = c("grey60", hhskblauw)) +
#   scale_y_continuous(limits = c(0,15), expand = c(0, 0)) +
#   labs(title = "Gemiddelde stikstofconcentratie",
#        subtitle = "per kwartaal",
#        caption = "De totale overschrijdingsfactor is de som van alle overschrijdingsfactoren in een monster",
#        y = "Overschrijdingsfactor",
#        x = "") +
#   theme(plot.subtitle = element_text(face = "plain"),
#         panel.grid.major.x = element_blank(),
#         axis.ticks.x = element_blank()) +
#   guides(fill = FALSE) +
#   NULL



data_n_plot <- data_n %>%
  group_by(jaarmaand, glas, tekst) %>%
  summarise(gem_waarde = mean(waarde)) %>%
  ungroup() %>%
  group_by(tekst) %>%
  filter(year(jaarmaand) > 2016)

n_plot <- data_n_plot %>%
  ggplot(aes(as.Date(jaarmaand), gem_waarde, colour = fct_rev(as.factor(glas)))) +
  annotate(x = c(q_begin, q_eind) -20, y = c(75, 75), geom = "area", fill = "grey85") +
  # geom_smooth(se = FALSE, linetype = "dashed", span = 0.8, size = 0.3) +
  geom_line(size = 1) +
  labs(x = "",
       y = "mg N/l",
       title = "Gemiddelde stikstofconcentratie per maand",
       subtitle = "Glastuinbouwgebied vs. overige gebieden",
       colour = "",
       caption = glue("De grijs gemarkeerde zone geeft {q_tekst} aan")) +
  scale_colour_manual(values = c("indianred3", hhskblauw), labels = c("Glastuinbouw", "Overig")) +
  hhskthema() +
  theme(plot.subtitle = element_text(size = 9)) +
  geom_text(aes(label = tekst), data = filter(data_n_plot, jaarmaand == max(jaarmaand)),
            x = q_eind + period(15, "days"), hjust = "left") +
  coord_cartesian(xlim = c(ymd(20161101, q_eind)), ylim = c(0,30), expand = FALSE, clip = "off") +
  guides(colour = FALSE) +
  theme(plot.margin = margin(5.5, 95, 5.5, 5.5, unit = "pt"))


data_p_plot <- data_p %>%
  group_by(jaarmaand, glas, tekst) %>%
  summarise(gem_waarde = mean(waarde)) %>%
  ungroup() %>%
  group_by(tekst) %>%
  filter(year(jaarmaand) > 2016)

p_plot  <-
  data_p_plot %>%
  ggplot(aes(as.Date(jaarmaand), gem_waarde, colour = fct_rev(as.factor(glas)))) +
  annotate(x = c(q_begin, q_eind) - 20, y = c(2.5, 2.5), geom = "area", fill = "grey85") +
  # geom_smooth(se = FALSE, linetype = "dashed", span = 0.8, size = 0.3) +
  geom_line(size = 1) +
  labs(x = "",
       y = "mg P/l",
       title = "Gemiddelde fosfaatconcentratie per maand",
       subtitle = "Glastuinbouwgebied vs. overige gebieden",
       colour = "",
       caption = glue("De grijs gemarkeerde zone geeft {q_tekst} aan")) +
  scale_colour_manual(values = c("indianred3", hhskblauw), labels = c("Glastuinbouw", "Overig")) +
  hhskthema() +
  theme(plot.subtitle = element_text(size = 9)) +
  geom_text(aes(label = tekst), data = filter(data_p_plot, jaarmaand == max(jaarmaand)),
            x = q_eind + period(15, "days"), hjust = "left") +
  coord_cartesian(xlim = ymd(20161101, q_eind), ylim = c(0, 2.5), expand = FALSE, clip = "off") +
  guides(colour = FALSE) +
  theme(plot.margin = margin(5.5, 95, 5.5, 5.5, unit = "pt"))




