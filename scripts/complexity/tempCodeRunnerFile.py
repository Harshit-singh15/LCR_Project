
    if most_freq_aa_percent < 50 and mutation_percent <= 50:
        return "CBR"

    elif most_freq_aa_percent >= 50 and mutation_percent <= 50:
        return "LCR"

    elif most_freq_aa_percent < 50 and mutation_percent > 50:
        return "HCR"