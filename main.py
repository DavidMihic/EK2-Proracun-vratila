from proracun import Vratilo, Lezaj

if __name__ == "__main__":
    # steps = [(diameter, length), ...] (mm)
    vratilo = Vratilo(
        "St 52-3",
        [(30, 20), (45, 45), (55, 110), (70, 24), (55, 110), (45, 41.5), (30, 35.5)],
        16,
        23,
    )
    vratilo.showDiagrams()

    vratilo.checkFlexuralCriticalRotationSpeed()
    print()
    vratilo.checkTorsionalCriticalRotationSpeed()

    print("\n")

    lezajA = Lezaj("NU 206 ECP", 30, vratilo.F_A, 0, 44000, 36500, 0)
    lezajA.checkBearing()

    print("\n")

    # First iteration
    lezajB = Lezaj("6406", 30, vratilo.F_B, vratilo.F_Ba, 43600, 23600, 12)
    lezajB.checkBearing()

    print()

    # Second iteration
    lezajB = Lezaj("6307", 35, vratilo.F_B, vratilo.F_Ba, 35100, 19000, 13)
    lezajB.checkBearing()
