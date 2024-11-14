## wget https://busco-data.ezlab.org/v5/data/info_mappings_all_busco_datasets_odb10.txt

awk -F '\t' ' BEGIN {print "mammalia_odb10 cetartiodactyla_odb10"} { if ($4=="cetartiodactyla_odb10") {C[$6]=$3} else { if ($4=="mammalia_odb10"&&$6 in C) {print $3,C[$6]} }} ' info_mappings_all_busco_datasets_odb10.txt
