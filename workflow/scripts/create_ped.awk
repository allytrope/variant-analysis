BEGIN { FS = "\t" }; {if ($5 == "M") SEX=1; else if ($5 == "F") SEX=2; else SEX=3}; {print $7,$1,$3,$4,SEX,"0"}
    