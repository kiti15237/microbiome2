

pie(c(10, 2, 28), labels = c("Constipation (10)", "Diarrhea (2)", "Normally formed (28)"), main="Control GI Dysfunction")
pie(c(11, 9, 18), labels = c("Constipation (11)", "Diarrhea (9)", "Normally formed (18)"), main="Autism GI Dysfunction")
x <- c(10, 2, 28)
y <- c(11, 9, 18)


aut_older <- sum(autMap$age_month_ok > cMap$age_month_ok)
control_older <- sum(autMap$age_month_ok < cMap$age_month_ok)
same_age <- sum(autMap$age_month_ok == cMap$age_month_ok)
                     
pie(c(aut_older, control_older, same_age ), 
    labels=c(paste("Aut older (", aut_older, ")"), paste("Control older (", control_older, ")"), paste("Same age (", same_age, ")")), 
    main="Older Sibling")
