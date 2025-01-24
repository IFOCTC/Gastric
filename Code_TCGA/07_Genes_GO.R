## Biological Process
angiogenesis <- c("COL18A1", "HIF1A", "SPARC", "NRP1", "CAV1",
                  "C3", "ITGB1", "MMP14", "CLIC4", "PGK1", "GRN")
blood_vessel_morphogenesis <- c("COL18A1", "HIF1A", "SPARC", "NRP1", "CAV1",
                                "C3", "ITGB1", "MMP14", "CLIC4", "PGK1", "GRN")
reg_cell_development <- c("DBN1", "GRN", "FBN1", "NRP1", "HIF1A",
                          "ITGB1", "GDI1", "TNFRSF1A", "MYADM")
positive_reg_cell_migration <- c("SPARC", "ITGB1", "MMP14", "NRP1", "HIF1A",
                                 "GNAI2", "GRN", "CAV1", "MYADM")
blood_vessel_development <- c("COL18A1", "HIF1A", "SPARC", "NRP1", "CAV1",
                              "C3", "ITGB1", "MMP14", "CLIC4", "PGK1", "GRN")
vasculature_development <- c("COL18A1", "HIF1A", "SPARC", "NRP1", "CAV1",
                             "C3", "ITGB1", "MMP14", "CLIC4", "PGK1", "GRN")
tube_morphogenesis <- c("COL18A1", "HIF1A", "SPARC", "NRP1", "CTSZ",
                        "CAV1", "C3", "ITGB1", "MMP14", "CLIC4", "PGK1", "GRN")
circulatory_system_development <- c("PDLIM1", "COL18A1", "HIF1A", "SPARC",
                                    "TNFRSF1A", "NRP1", "CAV1", "ECE1",
                                    "C3", "ITGB1", "MMP14", "CLIC4", "PGK1",
                                    "FBN1", "GRN")
epithelium_development <- c("COL18A1", "ITGB1", "GNAS", "NRP1", "HIF1A",
                            "CTSZ", "CAV1", "ECE1", "C3", "MMP14",
                            "TAF10", "CLIC4", "PGK1", "CD44", "TNFRSF1A",
                            "MYADM")
tube_development <- c("COL18A1", "HIF1A", "SPARC", "NRP1", "CTSZ",
                      "CAV1", "ECE1", "C3", "ITGB1", "MMP14", "LTBP3",
                      "CLIC4", "PGK1", "GRN")

## Hallmark
epithelial_mesenchymal_transition <- c("FBN1", "SPARC", "SERPINH1",
                                       "EMP3", "ITGB1", "CD44",
                                       "MMP14", "IGFBP4", "LGALS1")
coagulation <- c("C1S", "MMP14", "C1R", "C3", "SPARC", "FBN1", "FYN")
complement <- c("C1S", "C1R", "MMP14", "C3", "CTSD", "FYN", "GNAI2")

## Combine
bo_list       <- c(angiogenesis, blood_vessel_morphogenesis,
                      reg_cell_development, positive_reg_cell_migration,
                      blood_vessel_development, vasculature_development,
                      tube_morphogenesis, circulatory_system_development,
                      epithelium_development, tube_development)
bo_list_final <- unique(unique(bo_list))
hallmark_list <- c(epithelial_mesenchymal_transition, coagulation,
                      complement)
hallmark_list_final <- unique(unique(hallmark_list))

## Check
length(bo_list_final);bo_list_final
length(hallmark_list_final);hallmark_list_final

