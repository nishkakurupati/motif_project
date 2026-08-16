import motif_finder_utils_ppm
import evaluate_ppm_utils

string1 = "ACTGCA"
string2 = "TTTTGC"
string3 = "GACTCC"
sequences1 = ["AT",
              "GA",
              "GG",
              "CT"]
sequences2 = ["GT",
              "GT",
              "GG",
              "GT"]

motifs = [string1, string2, string3]

ppm1 = motif_finder_utils_ppm.build_ppm_from_brute_force(sequences1)
ppm2 = motif_finder_utils_ppm.build_ppm_from_brute_force(sequences2)

print(ppm2)

score = evaluate_ppm_utils.compare_ppm_get_score(ppm1, ppm2)
print(score)


