filter=50
for cell in GM12878 H1hESC IMR90; do
    rm "${cell}"/*-denoised.bin
    rm "${cell}"/*-filtered.bin
    for result in "${cell}"/*bin; do
        ~/phd/nQuire/nQuire view -c "${filter}" -o "${result/.bin/-filtered}" "${result}"
        ~/phd/nQuire/nQuire denoise -o "${result/.bin/-denoised}" "${result/.bin/-filtered}".bin
    done
done
~/phd/nQuire/nQuire lrdmodel -t 6 */*-chr*denoised.bin > nQuire-allResults.tsv
