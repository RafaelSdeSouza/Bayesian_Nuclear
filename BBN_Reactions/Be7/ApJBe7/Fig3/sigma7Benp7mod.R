sigma7Benp7mod <- function(ecm,
                           e0_1, ga_1, gb_1, ra_1, rb_1,
                           e0_2, ga_2, gb_2, ra_2, rb_2,
                           e0_3, ga_3, gb_3, ra_3, rb_3,
                           e0_4, ga_4, gb_4, ra_4, rb_4,
                           e0_5, ga_5, gb_5, ra_5, rb_5,
                           e0_6, ga_6, gb_6, ra_6, rb_6,
                           e0_7, ga_7, gb_7, ra_7, rb_7){
  
  SF1 <-  sigma7Benp(ecm, e0_1, ga_1, gb_1, ra_1, rb_1, jr = 2, la = 0, lb = 0)
  SF2 <-  sigma7Benp(ecm, e0_2, ga_2, gb_2, ra_2, rb_2, jr = 3, la = 1, lb = 1)
  SF3 <-  sigma7Benp(ecm, e0_3, ga_3, gb_3, ra_3, rb_3, jr = 3, la = 1, lb = 1)
  SF4 <-  sigma7Benp(ecm, e0_4, ga_4, gb_4, ra_4, rb_4, jr = 1, la = 0, lb = 0)
  SF5 <-  sigma7Benp(ecm, e0_5, ga_5, gb_5, ra_5, rb_5, jr = 4, la = 3, lb = 3)
  SF6 <-  sigma7Benp(ecm, e0_6, ga_6, gb_6, ra_6, rb_6, jr = 2, la = 1, lb = 1)
  SF7 <-   sigma7Benp(ecm, e0_7, ga_7, gb_7, ra_7, rb_7, jr = 0, la = 1, lb = 1)
  SF <- SF1 + SF2 + SF3 + SF4 + SF5 + SF6 + SF7
  return(SF = SF)
}