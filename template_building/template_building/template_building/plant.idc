﻿;IDA 4.80002 Form UTF-8
(DOCUMENT-HEADER :TYPE SCHEMA :PAGE-WIDTH 156 :PAGE-HEIGHT 112) 
(CONNECTION-LINE :AT ((449 428) (449 437) (498 437) (498 422) (528 422)) :FIRST-LINK ("EmeterCentChil" (0.464 1.0) OUTCONSUMLINK) :LAST-LINK ("ChilPw_Avg" (0.0 0.5) (|u| 2)) :DIR :RIGHT :ARROW (19 8 8)) 
(CONNECTION-LINE :AT ((413 429) (413 437) (498 437) (498 422) (528 422)) :FIRST-LINK ("EmeterLocalChil" (0.464 1) OUTCONSUMLINK) :LAST-LINK ("ChilPw_Avg" (0.0 0.5) (|u| 1)) :DIR :RIGHT :ARROW (19 8 8)) 
(CONNECTION-LINE :AT ((413 365) (413 353) (499 353) (499 371) (528 371)) :FIRST-LINK ("EmeterLocalBoil" (0.464 0.036) OUTCONSUMLINK) :LAST-LINK ("BoilPw_Avg" (0.0 0.536) (|u| 2)) :DIR :RIGHT :ARROW (19 8 8)) 
(CONNECTION-LINE :AT ((448 363) (448 353) (499 353) (499 371) (528 371)) :FIRST-LINK ("EmeterCentBoil" (0.429 0) OUTCONSUMLINK) :LAST-LINK ("BoilPw_Avg" (0.0 0.536) (|u| 1)) :DIR :RIGHT :ARROW (19 8 8)) 
(CONNECTION-LINE :AT ((552 422) (566 422)) :FIRST-LINK ("ChilPw_avg" (1 0.5) |uMean|) :LAST-LINK ("ChilPw_Snap" (0 0.5) (|u| 1)) :DIR :RIGHT :ARROW (19 8 8)) 
(CONNECTION-LINE :AT ((552 370) (566 370)) :FIRST-LINK ("BoilPw_Avg" (1 0.5) |uMean|) :LAST-LINK ("BoilPw_Snap" (0 0.5) (|u| 1)) :DIR :RIGHT :ARROW (19 8 8)) 
(EQUATION-FRAME :AT ((542 422)) :R (14.25 14.25) :ICON "lib:slidingaverage.ids" :SLOT ("ChilPw_Avg") :NAME "ChilPw_Avg" :DATA MODEL) 
(EQUATION-FRAME :AT ((586 422)) :R (14.25 14.25) :ICON "lib:snapminmax.ids" :SLOT ("ChilPw_Snap") :NAME "ChilPw_Snap" :DATA MODEL) 
(EQUATION-FRAME :AT ((542 370)) :R (14.25 14.25) :ICON "lib:slidingaverage.ids" :SLOT ("BoilPw_Avg") :NAME "BoilPw_Avg" :DATA MODEL) 
(EQUATION-FRAME :AT ((586 370)) :R (14.25 14.25) :ICON "lib:snapminmax.ids" :SLOT ("BoilPw_Snap") :NAME "BoilPw_Snap" :DATA MODEL) 
(TEXT-OBJECT :VALUE "DHW control" :AT ((239 49) (285 87)) :STYLE NOTE) 
(CONNECTION-LINE :AT ((234 44) (234 176) (251 176)) :LINE-STYLE :DOT :FIRST-LINK (:SELF (0.369 0.0) "DHW_ctl") :LAST-LINK ("boil" (0 0.214) DOMWATCONTROL) :DIR :RIGHT :ARROW (19 8 8)) 
(CONNECTION-LINE :AT ((406 219) (406 188) (449 188)) :FIRST-LINK ("TAmbRef1" (0.5 0.107) LINK) :LAST-LINK ("chil" (0 0.429) TCONDENS) :DIR :RIGHT :ARROW (19 8 8)) 
(LINK-FRAME :AT ((376 216) (436 244)) :STYLE HYPERLINK :SLOT ("TAmbRef1") :NAME "TAmbRef1" :SIDE :TOP :LABEL "TAir2" :ICON CLIMATE-PIXMAP :DATA REFERENCE) 
(CONNECTION-LINE :AT ((208 241) (332 241) (332 173) (308 173)) :LINE-STYLE :DOT :FIRST-LINK ("Plantctr_Sched" (1.0 0.531) OUTSIGNALLINK) :LAST-LINK ("boil" (1.0 0.161) PUMPCONTROL) :DIR :RIGHT :ARROW (19 8 8)) 
(TEXT-OBJECT :VALUE (GERMAN "Wärme- und Kälteerzeuger mit sehr hoher Kapazität. Heißwasservlorlauftemperatur ist abhängig von der Außentemperatur, Kaltwasser ist konstant.
Wärme- und Kälteerzeuger enthalten weitere Parameter." ENGLISH "Plant model with (by default) very large capacity. Supply hot water setpoint is a function of outside air temp. Chilled water temperatures to zones and AHU are constant.
Open boiler and chiller to set parameters." FINNISH "Rajoittamaton tehon (oletus) lämmöntuoton malli. Lämmityksen menoveden lämpötilan asetusarvo riippuu ulkoilman lämpötilasta. Jäähdytetyn veden lämpötilat vyöhykkeisiin ja IV-koneeseen ovat vakioita.
Avaa primäärijärjestelmä syöttääksesi parametrit." FRENCH "Modèle de centrale à très grande capacité (par défaut). La valeur de consigne de la température d’eau chaude est une fonction de la température extérieure. Températures d’eau froide vers les zones et la CTA sont constantes.
Ouvrir la chaudière et la machine frigorifique pour fixer les paramètres." SWEDISH "Energicentral med stor kapacitet. Börvärdet för framledningstemperaturen beror av utomhustemperaturen. Köldbärarens temperatur till zoner respektive LB är konstant.
Öppna värmepanna respektive kylmaskin för att ställa in egenskaper." SPANISH "Modelo de planta central con (por defecto) una capacidad muy grande. Punto de consigna del agua caliente de suministro está en función de la temp. del aire exterior. Las temperaturas del agua refrigerada hacia las zonas y las UTAs son constantes.
Abrir la caldera y la enfriadora para establecer los parámetros.") :AT ((20 360) (320 428)) :STYLE NOTE :LINE-GROUND #S(RGB RED 255 GREEN 255 BLUE 233) :OUTLINED-P T :MARKUP HTML :PADDING 5) 
(CONNECTION-LINE :AT ((164 169) (178 169) (178 138)) :FIRST-LINK ("OffSet" (1.0 0.531) OUTSIGNALLINK) :LAST-LINK ("ADDER" (0.25 1) (INSIGNALLINK 1)) :DIR :RIGHT :ARROW (19 8 8)) 
(TEXT-OBJECT :VALUE (ENGLISH "Standard Plant" FINNISH "Oletus lämmön- ja jäädytyksentuotto" GERMAN "Wärme- und Kälteerzeuger" FRENCH "Centrale standard" SWEDISH "Standardanläggning" SPANISH "Planta central estándar") :TEXT-COLOR :DEFAULT :FONT (:SWISS :ARIAL 16 1) :AT ((16 12) (379 32)) :STYLE SECTION) 
(CONNECTION-LINE :AT ((457 103) (493 103) (493 164)) :LINE-STYLE :DOT :FIRST-LINK ("Sched_Plant" (1.0 0.531) OUTSIGNALLINK) :LAST-LINK ("chil" 3.268 PUMPCONTROL) :DIR :RIGHT :ARROW (19 8 8)) 
(CONNECTION-LINE :AT ((65 153) (65 116) (81 116)) :FIRST-LINK ("TAmbRef" (0.483 0.179) LINK) :LAST-LINK ("Setp_Boil" (0 0.563) INSIGNALLINK) :DIR :RIGHT :ARROW (19 8 8)) 
(CONNECTION-LINE :AT ((134 116) (148 116) (166 116)) :FIRST-LINK ("Setp_Boil" (1.0 0.563) OUTSIGNALLINK) :LAST-LINK ("ADDER" (0 0.5) (INSIGNALLINK 2)) :DIR :RIGHT :ARROW (19 8 8)) 
(FRAME-BOX :VALUE (ENGLISH "Setpoint for supply hot water" GERMAN "Sollwert für Heißwasser" FINNISH "Lämmityksen menoveden lämpötilan asetusarvo" FRENCH "Valeur de consigne pour l’eau chaude fournie" SWEDISH "Börvärde för värmesystemets framledningstemperatur" SPANISH "Valores de consigna para el agua caliente de suministro") :AT ((36 60) (224 208)) :STYLE SECTION) 
(CONNECTION-LINE :AT ((213 116) (320 116) (320 185) (305 185)) :FIRST-LINK ("ADDER" (1 0.5) OUTSIGNALLINK) :LAST-LINK ("boil" (1 0.375) (TEMPSETP 1)) :DIR :RIGHT :ARROW (19 8 8)) 
(CONNECTION-LINE :AT ((213 116) (328 116) (328 185) (305 185)) :FIRST-LINK ("ADDER" (1 0.5) OUTSIGNALLINK) :LAST-LINK ("boil" (1 0.375) (TEMPSETP 2)) :DIR :RIGHT :ARROW (19 8 8)) 
(LINK-FRAME :AT ((36 148) (96 176)) :STYLE HYPERLINK :SLOT (:BUILDING CLIMATE TAIR2) :NAME "TAmbRef" :SIDE :TOP :LABEL "TAir2" :ICON CLIMATE-PIXMAP :DATA REFERENCE) 
(CONNECTION-LINE :AT ((308 202) (351 202) (351 287) (595 287)) :LINE-COLOR (:CALL PMT-COLOR [@ 1] [@ 2]) :LINE-STYLE 3 :FIRST-LINK ("boil" (1.0 0.679) (INLET 2)) :LAST-LINK (:SELF 2.191 "Zone_rtn_hot")) 
(CONNECTION-LINE :AT ((310 197) (371 197) (371 250) (595 250)) :LINE-COLOR (:CALL PMT-COLOR [@ 1] [@ 2]) :LINE-STYLE 3 :FIRST-LINK ("boil" (1 0.589) (OUTLET 2)) :LAST-LINK (:SELF 2.312 "Zone_sup_hot")) 
(CONNECTION-LINE :AT ((310 197) (348 197) (348 44)) :LINE-COLOR (:CALL PMT-COLOR [@ 1] [@ 2]) :LINE-STYLE 3 :FIRST-LINK ("boil" (1 0.589) (OUTLET 1)) :LAST-LINK (:SELF (0.573 0.0) "AHU_sup_hot")) 
(CONNECTION-LINE :AT ((308 202) (337 202) (337 44)) :LINE-COLOR (:CALL PMT-COLOR [@ 1] [@ 2]) :LINE-STYLE 3 :FIRST-LINK ("boil" (1.0 0.679) (INLET 1)) :LAST-LINK (:SELF 3.446 "AHU_rtn_hot")) 
(CONNECTION-LINE :AT ((508 176) (532 176) (532 44)) :LINE-COLOR (:CALL PMT-COLOR [@ 1] [@ 2]) :LINE-STYLE 3 :FIRST-LINK ("chil" (1.0 0.214) OUTLET1) :LAST-LINK (:SELF (0.891 0.0) "AHU_sup_cold")) 
(CONNECTION-LINE :AT ((508 186) (545 186) (545 44)) :LINE-COLOR (:CALL PMT-COLOR [@ 1] [@ 2]) :LINE-STYLE 3 :FIRST-LINK ("chil" (1.0 0.393) INLET1) :LAST-LINK (:SELF 3.09 "AHU_rtn_cold")) 
(CONNECTION-LINE :AT ((508 204) (508 205) (549 205) (595 205)) :LINE-COLOR (:CALL PMT-COLOR [@ 1] [@ 2]) :LINE-STYLE 3 :FIRST-LINK ("chil" (1.0 0.714) OUTLET2) :LAST-LINK (:SELF 2.463 "Zone_sup_cold")) 
(CONNECTION-LINE :AT ((509 212) (558 212) (596 212)) :LINE-COLOR (:CALL PMT-COLOR [@ 1] [@ 2]) :LINE-STYLE 3 :FIRST-LINK ("chil" (1 0.857) INLET2) :LAST-LINK (:SELF (0.995 0.56) "Zone_rtn_cold")) 
(EQUATION-FRAME :AT ((374 414)) :R (14 14) :ICON "lib:emeter.ids" :SLOT ("EmeterWater") :NAME "EmeterWater" :DATA :CEO :D "Calculates the total power consumption and the cost of energy per hour for a number of monitored input power signals. A signal may be multiplied by a factor (mult), signifying multiple consumers.") 
(EQUATION-FRAME :AT ((374 378)) :R (14 14) :ICON "lib:emeter.ids" :SLOT ("EmeterPump") :NAME "EmeterPump" :DATA :CEO :D "Calculates the total power consumption and the cost of energy per hour for a number of monitored input power signals. A signal may be multiplied by a factor (mult), signifying multiple consumers.") 
(EQUATION-FRAME :AT ((450 414)) :R (14 14) :ICON "lib:emeter.ids" :SLOT ("EmeterCentChil") :NAME "EmeterCentChil" :DATA :CEO :D "Calculates the total power consumption and the cost of energy per hour for a number of monitored input power signals. A signal may be multiplied by a factor (mult), signifying multiple consumers.") 
(EQUATION-FRAME :AT ((450 378)) :R (14 14) :ICON "lib:emeter.ids" :SLOT ("EmeterCentBoil") :NAME "EmeterCentBoil" :DATA :CEO :D "Calculates the total power consumption and the cost of energy per hour for a number of monitored input power signals. A signal may be multiplied by a factor (mult), signifying multiple consumers.") 
(EQUATION-FRAME :AT ((414 414)) :R (14 14) :ICON "lib:emeter.ids" :SLOT ("EmeterLocalChil") :NAME "EmeterLocalChil" :DATA :CEO :D "Calculates the total power consumption and the cost of energy per hour for a number of monitored input power signals. A signal may be multiplied by a factor (mult), signifying multiple consumers.") 
(EQUATION-FRAME :AT ((414 378)) :R (14 14) :ICON "lib:emeter.ids" :SLOT ("EmeterLocalBoil") :NAME "EmeterLocalBoil" :DATA :CEO :D "Calculates the total power consumption and the cost of energy per hour for a number of monitored input power signals. A signal may be multiplied by a factor (mult), signifying multiple consumers.") 
(EQUATION-FRAME :AT ((192 116)) :R (20 20) :ICON "lib:adder.ids" :SLOT ("ADDER") :NAME "ADDER" :DATA :EO :D "Adds n input signals. Sends sum to multiple output links.") 
(EQUATION-FRAME :AT ((144 168)) :R (24 16) :ICON "lib:schedule.ids" :SLOT ("OffSet") :TITLE "Setback schedule" :NAME "OffSet" :DATA SCHEDULE :HELP-STRING (ENGLISH "Schedule for boiler setpoint offset" FINNISH "Aikataulu lämmöntuoton asetusarvon muutokselle" GERMAN "Zeitplan für Differenz in Heizkesselsollwert" FRENCH "Horaire pour le décalage de la valeur de consigne de la chaudière" SWEDISH "Tidsschema för parallellförflyttning av värmepannans börvärde" SPANISH "Horario para el desfase del punto de consigna de la caldera")) 
(EQUATION-FRAME :AT ((184 240)) :R (24 16) :ICON "lib:schedule.ids" :SLOT ("Plantctr_Sched") :TITLE "Boiler operation" :NAME "Plantctr_Sched" :DATA SCHEDULE :HELP-STRING (ENGLISH "Schedule for boiler operation" FINNISH "Aikataulu lämmöntuotolle" GERMAN "Zeitplan für Betrieb Heizkesselsollwert" FRENCH "Horaire d’exploitation de la chaudière" SWEDISH "Tidsschema för värmepannans drift" SPANISH "Horario de funcionamiento de la caldera")) 
(EQUATION-FRAME :AT ((437 102)) :R (24 16) :ICON "lib:schedule.ids" :SLOT ("Sched_Plant") :TITLE "Chiller operation" :NAME "Sched_Plant" :DATA SCHEDULE :HELP-STRING (ENGLISH "Algorithmic schedule object for IDA ICE; leap years and DST are handled" FINNISH "Algoritminen objekti aikataululle IDA Indoor Climate and Energy-sovelluksessa. Karkausvuosi ja kesäaika otetaan huomioon." GERMAN "Dies ist ein algorithmischer Zeitplan für IDA ICE. Schaltjahre und Sommerzeit werden berücksichtigt." NORWEGIAN "Dette er algoritmeobjekt for urstyring i IDA Klima og Energi. Det tas hensyn til skuddår og sommertid" FRENCH "Objet horaire algorithmique pour IDA ICE. Les années bissextiles et l’heure d’été y sont gérées." SWEDISH "Algoritmiskt tidsschema; skottår och sommartid hanteras" SPANISH "Objeto horario algorítmico para IDA ICE ; anos bisiestos y el horario de verano son gestionados")) 
(EQUATION-FRAME :AT ((110 114)) :R (24 16) :ICON "lib:PLINSEGM.ids" :SLOT ("Setp_Boil") :NAME "Setp_Boil" :DATA :EO :HELP-STRING (ENGLISH "Supply temperature setpoint" FINNISH "IV-kon. menoveden lämpötila" GERMAN "Sollwert für Zulufttemperatur" NORWEGIAN "Børverdi for turtemperatur" FRENCH "Valeur de consigne de la température fournie" SWEDISH "Börvärde för framledningstemperatur" SPANISH "Punto de consigna de la temperatura de suministro")) 
(EQUATION-FRAME :AT ((480 192)) :R (28 28) :ICON "lib:simchil.ids" :SLOT ("chil") :NAME "chil" :DATA :EO :HELP-STRING (ENGLISH "A simplified chiller and water pump model with ideal control." FINNISH "Yksinkertaistettu malli jäähdytyskoneelle ja vesipiiri ideaalisäädöllä" GERMAN "Vereinfachtes Modell für Kältemaschine und Pumpenkreis mit idealer Regelung." NORWEGIAN "Forenklet modell av kjølemaskin og pumpekrets med ideell regulering" FRENCH "Un modèle simplifié de machine frigorifique et pompe à eau avec réglage idéal." SWEDISH "Förenklad kylmaskin och pump med idealiserad styrning" SPANISH "Un modelo simplificado de enfriadora y bomba de agua con control ideal.")) 
(EQUATION-FRAME :AT ((280 192)) :R (28 28) :ICON "lib:simboil.ids" :SLOT ("boil") :NAME "boil" :DATA :EO :HELP-STRING (ENGLISH "Simple boiler with pump circuit" FINNISH "Yksinkertaistettu malli lämmöntuotolle ja vesipiirille" GERMAN "Vereinfachtes Modell für Heizkessel und Heizkreis" NORWEGIAN "Enkel modell for kjel med pumpekrets" FRENCH "Chaudière simple avec circuit de pompe" SWEDISH "Enkel värmepanna med pump" SPANISH "Caldera simple con circuito de bombeo")) 
(SELF-FRAME :AT ((310 194)) :R (289 150.0) :SLOT (:SELF) :DATA MACRO-OBJECT) 
;[end of template_building\plant.idc]