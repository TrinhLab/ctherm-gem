rem=$(cat imbalances_other.csv | sed 1d | cut -d ',' -f1)
for i in "${rem[@]}"
do
	grep "$i" imbalances_remaining.csv
done
