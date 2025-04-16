import numpy as np
from math import sqrt, cos, tan, radians as rad
import matplotlib.pyplot as plt


def linearInterpolate(x, xRanges, yRanges):
    rightIdx = [i for i in range(len(xRanges)) if xRanges[i] >= x][0]

    if rightIdx == 0:
        return yRanges[0]

    xLeft = xRanges[rightIdx - 1]
    xRight = xRanges[rightIdx]
    yLeft = yRanges[rightIdx - 1]
    yRight = yRanges[rightIdx]

    return yLeft + (x - xLeft) / (xRight - xLeft) * (yRight - yLeft)


class Vratilo:
    middlePartWidth = 24  # mm

    sigmaF_max = 50  # N/mm^2
    sigma_fDN = 240  # N/mm^2
    tau_fDI = 190  # N/mm^2
    Rm = 500  # N/mm^2

    ALPHA = 20  # deg
    ALPHA_N = 20  # deg
    BETA = 18  # deg

    resolution = 0.5  # mm

    def __init__(self, steps, leftBearingWidth, rightBearingWidth):
        self.T = 420000  # Nmm
        self.J2 = 0.0500  # kgm^2 = Nms^2
        self.J3 = 1.050  # kgm^2 = Nms^2
        self.l = 350  # mm
        self.G_z2 = 210  # N
        self.G_z3 = 80  # N
        self.b2 = 110  # mm
        self.b3 = 110  # mm
        self.r2 = 165  # mm
        self.r3 = 57.8  # mm

        self.steps = steps
        self.L = sum(pair[1] for pair in steps)
        self.leftBearingWidth = leftBearingWidth
        self.rightBearingWidth = rightBearingWidth

        self.l3 = (self.l - self.b2 - self.middlePartWidth) / 2
        self.l6 = (self.l + self.b3 + self.middlePartWidth) / 2

        self.criticalSections = [
            leftBearingWidth / 2,
            steps[1][1] + leftBearingWidth / 2,
            self.l3,
            self.l3 + self.b2 / 2,
            self.l6 - self.b3 / 2,
            self.l6,
            self.l - steps[-2][1] - rightBearingWidth / 2,
            self.l - rightBearingWidth / 2,
        ]

        self.F_t2 = self.T / self.r2
        self.F_r2 = self.F_t2 * tan(rad(self.ALPHA))
        self.T2 = self.T

        self.F_t3 = self.T / self.r3
        self.F_r3 = self.F_t3 * tan(rad(self.ALPHA_N)) / cos(rad(self.BETA))
        self.F_a3 = self.F_t3 * tan(rad(self.BETA))
        self.T3 = self.T

        self.F_Ah = (
            self.F_r3 * (self.l - self.l6)
            + self.F_a3 * self.r3
            - self.F_r2 * (self.l - self.l3)
        ) / self.l
        self.F_Av = (
            (self.F_t2 + self.G_z2) * (self.l - self.l3)
            + (self.F_t3 + self.G_z3) * (self.l - self.l6)
        ) / self.l
        self.F_A = sqrt(self.F_Ah**2 + self.F_Av**2)

        self.F_Bh = -self.F_Ah - self.F_r2 + self.F_r3
        self.F_Bv = self.F_t2 + self.G_z2 + self.F_t3 + self.G_z3 - self.F_Av
        self.F_Ba = self.F_a3
        self.F_B = sqrt(self.F_Bh**2 + self.F_Bv**2)

        x = list(np.arange(0, self.l + self.resolution, self.resolution))
        l3_index = [i for i in range(len(x)) if x[i] == self.l3][0]
        x = x[: l3_index + 1] + [self.l3] + x[l3_index + 1 :]
        l6_index = [i for i in range(len(x)) if x[i] == self.l6][0]
        self.x = x[: l6_index + 1] + [self.l6] + x[l6_index + 1 :]
        self.x1 = self.x[: l3_index + 1]
        self.x2 = self.x[l3_index + 1 : l6_index + 1]
        self.x3 = self.x[l6_index + 1 :]

        self.Nx1 = lambda x: 0
        self.Nx2 = lambda x: 0
        self.Nx3 = lambda x: self.F_a3

        self.Qy1 = lambda x: self.F_Ah
        self.Qy2 = lambda x: self.F_Ah + self.F_r2
        self.Qy3 = lambda x: -self.F_Bh

        self.Qz1 = lambda x: self.F_Av
        self.Qz2 = lambda x: self.F_Av - self.F_t2 - self.G_z2
        self.Qz3 = lambda x: -self.F_Bv

        self.Tx1 = lambda x: 0
        self.Tx2 = lambda x: -self.T
        self.Tx3 = lambda x: 0

        self.My1 = lambda x: self.F_Av * x
        self.My2 = lambda x: (self.F_Av - self.F_t2 - self.G_z2) * (
            x - self.l3
        ) + self.My1(self.l3)
        self.My3 = lambda x: -self.F_Bv * (x - self.l6) + self.My2(self.l6)

        self.Mz1 = lambda x: -self.F_Ah * x
        self.Mz2 = lambda x: -(self.F_Ah + self.F_r2) * (x - self.l3) + self.Mz1(
            self.l3
        )
        self.Mz3 = (
            lambda x: self.F_Bh * (x - self.l6)
            + self.F_a3 * self.r3
            + self.Mz2(self.l6)
        )

        self.alpha0 = self.sigma_fDN / (sqrt(3) * self.tau_fDI)
        self.k = (10 / self.sigmaF_max) ** (1 / 3)

    def showDiagrams(self, darkMode=True):
        self.plotForcesMoments(darkMode)
        self.plotShaft(darkMode)
        plt.show()

    def plotShaft(self, darkMode=True, idealShaftLineColor="red"):
        plt.style.use("dark_background" if darkMode else "default")
        plt.figure("Idealni oblik vratila", figsize=(13, 8))
        plt.tight_layout(pad=2, h_pad=5)
        plt.suptitle("Idealni oblik vratila")
        plt.grid(True, alpha=0.3)

        plt.axhline(0, color="white" if darkMode else "black", ls="-.")

        radii = self.getRadius(self.x)

        # Plot ideal shaft
        plt.plot(self.x, radii, idealShaftLineColor)
        plt.plot(self.x, [-r for r in radii], idealShaftLineColor)

        extraLinesColor = "white" if darkMode else "black"

        # Plot steps
        curX = self.leftBearingWidth / 2 - self.steps[0][1]
        prevRadius = 0
        for diameter, distance in self.steps + [(0, 0)]:
            radius = diameter / 2
            plt.plot([curX, curX + distance], [radius, radius], color=extraLinesColor)
            plt.plot([curX, curX + distance], [-radius, -radius], color=extraLinesColor)
            plt.plot([curX, curX], [prevRadius, radius], color=extraLinesColor)
            plt.plot([curX, curX], [-prevRadius, -radius], color=extraLinesColor)
            curX += distance
            prevRadius = radius

        # Indicate critical sections
        maxRadius = max(self.steps)[0] / 2
        for i, distance in enumerate(self.criticalSections):
            plt.plot(
                [distance, distance],
                [-maxRadius * 1.1, maxRadius * 1.1],
                color=extraLinesColor,
                ls="--",
            )
            plt.annotate(
                str(i + 1),
                xy=(distance, maxRadius * 1.15),
                xytext=(distance - 2.5, maxRadius * 1.15),
                color="white",
                fontsize=15,
                va="center",
            )
            plt.annotate(
                str(i + 1),
                xy=(distance, -maxRadius * 1.175),
                xytext=(distance - 2.5, -maxRadius * 1.175),
                color="white",
                fontsize=15,
                va="center",
            )

    def getDiameter(self, x):
        if isinstance(x, (int, float)):
            if x < self.l3:
                M = sqrt(self.Mz1(x) ** 2 + self.My1(x) ** 2)
            elif x <= self.l6:
                Mf = sqrt(self.Mz2(x) ** 2 + self.My2(x) ** 2)
                M = sqrt(Mf**2 + 0.75 * (self.alpha0 * self.Tx2(x)) ** 2)
            else:
                M = sqrt(self.Mz3(x) ** 2 + self.My3(x) ** 2)

            return self.k * M ** (1 / 3)

        d = [0] * len(x)
        for i in range(len(x)):
            if (x[i] == self.l3 and x[i] == x[i + 1]) or (x[i] < self.l3):
                M = sqrt(self.Mz1(x[i]) ** 2 + self.My1(x[i]) ** 2)
            elif (x[i] == self.l6 and x[i] == x[i + 1]) or (x[i] < self.l6):
                Mf = sqrt(self.Mz2(x[i]) ** 2 + self.My2(x[i]) ** 2)
                M = sqrt(Mf**2 + 0.75 * (self.alpha0 * self.Tx2(x[i])) ** 2)
            else:
                M = sqrt(self.Mz3(x[i]) ** 2 + self.My3(x[i]) ** 2)
            d[i] = self.k * M ** (1 / 3)

        return d

    def getRadius(self, x):
        if isinstance(x, (int, float)):
            return self.getDiameter(x) / 2

        return [d / 2 for d in self.getDiameter(x)]

    def plotForcesMoments(self, darkMode=True):
        plt.style.use("dark_background" if darkMode else "default")
        fig, axes = plt.subplots(2, 3, figsize=(15, 10))
        fig.suptitle("Sile i momenti")
        fig.canvas.manager.set_window_title("Sile i momenti")

        graphs = [
            [[self.Nx1, self.Nx2, self.Nx3], "Nx", "Sila (kN)"],
            [[self.Qy1, self.Qy2, self.Qy3], "Qy", "Sila (kN)"],
            [[self.Qz1, self.Qz2, self.Qz3], "Qz", "Sila (kN)"],
            [[self.Tx1, self.Tx2, self.Tx3], "Tx", "Moment (Nm)"],
            [[self.My1, self.My2, self.My3], "My", "Moment (Nm)"],
            [[self.Mz1, self.Mz2, self.Mz3], "Mz", "Moment (Nm)"],
        ]

        for i, (y, title, label) in enumerate(graphs):
            ax = axes[i // 3, i % 3]
            self.plotDiagram(y, ax, darkMode)
            ax.grid(True, alpha=0.3)
            ax.set_title(title)
            ax.set_xlabel("Udaljenost (mm)")
            ax.set_ylabel(label)

        plt.tight_layout(rect=[0, 0, 1, 0.95])

    def plotDiagram(self, y, ax, darkMode, lineColor="aquamarine"):
        extraLineColor = "white" if darkMode else "black"
        ax.axhline(0, color=extraLineColor)
        ax.axvline(0, linestyle="--", color=extraLineColor)
        ax.axvline(self.l3, linestyle="--", color=extraLineColor)
        ax.axvline(self.l6, linestyle="--", color=extraLineColor)
        ax.axvline(self.l, linestyle="--", color=extraLineColor)

        xData = [self.x1, self.x2, self.x3]
        for i in range(len(xData)):
            curX = xData[i]
            curY = y[i]
            yValues = [curY(x) / 1000 for x in curX]  # N -> kN, Nmm -> Nm
            ax.plot(curX, yValues, color=lineColor, linewidth=2)
            ax.fill_between(
                curX,
                yValues,
                0,
                where=[val >= 0 for val in yValues],
                color="green",
                alpha=0.3,
            )
            ax.fill_between(
                curX,
                yValues,
                0,
                where=[val < 0 for val in yValues],
                color="red",
                alpha=0.3,
            )

    def getReducedMoment(self, x):
        if x < self.l3:
            return sqrt(self.Mz1(x) ** 2 + self.My1(x) ** 2)
        elif x <= self.l6:
            return sqrt(self.Mz2(x) ** 2 + self.My2(x) ** 2)
        return sqrt(self.Mz3(x) ** 2 + self.My3(x) ** 2)


class Lezaj:
    n_m = 320  # rpm
    L_h = 8000  # hrs
    helperRatioRanges = [0.172, 0.345, 0.689, 1.03, 1.38, 2.07, 3.45, 5.17, 6.89]
    eRanges = [0.19, 0.22, 0.26, 0.28, 0.30, 0.34, 0.38, 0.42, 0.44]
    yRanges = [2.30, 1.99, 1.71, 1.55, 1.45, 1.31, 1.15, 1.04, 1.00]

    def __init__(
        self,
        designation,
        innerDiameter,
        radialForce,
        axialForce,
        C,
        C0,
        f0,
        isLeftBearing,
    ):
        self.designation = designation
        self.innerDiameter = innerDiameter
        self.radialForce = radialForce
        self.axialForce = axialForce
        self.forceRatio = axialForce / radialForce
        self.C = C
        self.C0 = C0
        self.f0 = f0
        self.isLeftBearing = isLeftBearing

        if isLeftBearing:
            self.epsilon = 10 / 3
        else:
            self.epsilon = 3

    def checkBearing(self):
        self.helperRatio = self.f0 * self.axialForce / self.C0
        self.e = linearInterpolate(
            self.helperRatio, self.helperRatioRanges, self.eRanges
        )

        if self.isLeftBearing or self.forceRatio <= self.e:
            self.X, self.Y = 1, 0
        else:
            self.X = 0.56
            self.Y = linearInterpolate(self.e, self.eRanges, self.yRanges)

        radialLoadEquivalent = self.X * self.radialForce + self.Y * self.axialForce

        self.C1 = radialLoadEquivalent * (60 * self.n_m * self.L_h / 1e6) ** (
            1 / self.epsilon
        )
        self.Lh = (
            1e6 / (60 * self.n_m) * (self.C / radialLoadEquivalent) ** self.epsilon
        )

        print("Ležaj:", self.designation)
        print("C1 (N):", self.C1)
        print("Vrijeme (h):", self.Lh)

        print(f"{'' if self.C1 <= self.C else 'NE '}ZADOVOLJAVA!")
