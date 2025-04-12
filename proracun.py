import numpy as np
from math import sqrt, cos, tan, radians as rad
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon


class Vratilo:
    middlePartWidth = 24  # mm

    sigmaF_max = 50  # N/mm^2
    sigma_fDN = 240  # N/mm^2
    tau_fDI = 190  # N/mm^2

    ALPHA = 20  # deg
    ALPHA_N = 20  # deg
    BETA = 18  # deg

    resolution = 0.5  # mm

    def __init__(self, debug=False):
        if debug:  # Sample values
            self.T = 610000  # Nmm
            self.J2 = 0.0500  # kgm^2 = Nms^2
            self.J3 = 0.8750  # kgm^2 = Nms^2
            self.l = 360  # mm
            self.G_z2 = 240  # N
            self.G_z3 = 110  # N
            self.b2 = 120  # mm
            self.b3 = 120  # mm
            self.r2 = 180  # mm
            self.r3 = 63.1  # mm
        else:
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

        self.l3 = (self.l - self.b2 - self.middlePartWidth) / 2
        self.l6 = (self.l + self.b3 + self.middlePartWidth) / 2

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

        x = list(np.arange(0, self.l, self.resolution))
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

    def plotIdealShaft(self):
        plt.figure("Idealni oblik vratila", figsize=(10, 7))
        plt.suptitle("Idealni oblik vratila")
        plt.grid(True)

        plt.axhline(0, color="black")
        plt.axvline(0, linestyle="--", color="black")
        plt.axvline(self.l3, linestyle="--", color="black")
        plt.axvline(self.l6, linestyle="--", color="black")
        plt.axvline(self.l, linestyle="--", color="black")

        diameters = self.getDiameter(self.x)
        plt.plot(self.x, diameters, "b")

        vertices = [(self.x[0], 0)] + list(zip(self.x, diameters)) + [(self.x[-1], 0)]
        polygon = Polygon(
            vertices,
            closed=True,
            facecolor="none",
            edgecolor="b",
            hatch="//",
        )
        plt.gca().add_patch(polygon)

        plt.xlabel("Udaljenost (mm)")
        plt.ylabel("Promjer vratila (mm)")

        plt.show()

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

    def plotForcesMoments(self):
        fig, axes = plt.subplots(2, 3, figsize=(15, 10))
        fig.suptitle("Sile i momenti")
        fig.canvas.manager.set_window_title("Sile i momenti")

        graphs = [
            [[self.Nx1, self.Nx2, self.Nx3], "y", "Nx", "Sila (N)"],
            [[self.Qy1, self.Qy2, self.Qy3], "r", "Qy", "Sila (N)"],
            [[self.Qz1, self.Qz2, self.Qz3], "c", "Qz", "Sila (N)"],
            [[self.Tx1, self.Tx2, self.Tx3], "b", "Tx", "Moment (Nmm)"],
            [[self.My1, self.My2, self.My3], "m", "My", "Moment (Nmm)"],
            [[self.Mz1, self.Mz2, self.Mz3], "g", "Mz", "Moment (Nmm)"],
        ]

        for i, (y, color, title, label) in enumerate(graphs):
            ax = axes[i // 3, i % 3]
            self.plotDiagram(y, color, ax)
            ax.grid(True)
            ax.set_title(title)
            ax.set_xlabel("Udaljenost (mm)")
            ax.set_ylabel(label)

        plt.tight_layout(rect=[0, 0, 1, 0.95])
        plt.show()

    def plotDiagram(self, y, color, ax):
        ax.axhline(0, color="black")
        ax.axvline(0, linestyle="--", color="black")
        ax.axvline(self.l3, linestyle="--", color="black")
        ax.axvline(self.l6, linestyle="--", color="black")
        ax.axvline(self.l, linestyle="--", color="black")

        xData = [self.x1, self.x2, self.x3]
        for i in range(len(xData)):
            curX = xData[i]
            curY = y[i]
            ax.plot(curX, list(map(curY, curX)), color=color, linewidth=2)
            ax.fill_between(curX, list(map(curY, curX)), 0, color=color, alpha=0.3)


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
        self.rightIdx = [
            i
            for i in range(len(self.helperRatioRanges))
            if self.helperRatioRanges[i] >= self.helperRatio
        ][0]

        helperRatioLeft = self.helperRatioRanges[self.rightIdx - 1]
        helperRatioRight = self.helperRatioRanges[self.rightIdx]
        eLeft = self.eRanges[self.rightIdx - 1]
        eRight = self.eRanges[self.rightIdx]

        self.e = eLeft + (self.helperRatio - helperRatioLeft) / (
            helperRatioRight - helperRatioLeft
        ) * (eRight - eLeft)

        if self.isLeftBearing or self.forceRatio <= self.e:
            self.X, self.Y = 1, 0
        else:
            yLeft = self.yRanges[self.rightIdx - 1]
            yRight = self.yRanges[self.rightIdx]
            self.X = 0.56
            self.Y = yLeft + (self.e - eLeft) / (eLeft - eRight) * (yRight - yLeft)

        radialLoadEquivalent = self.X * self.radialForce + self.Y * self.axialForce

        self.C1 = radialLoadEquivalent * (60 * self.n_m * self.L_h / 1e6) ** (
            1 / self.epsilon
        )
        self.Lh = (
            1e6 / (60 * self.n_m) * (self.C / radialLoadEquivalent) ** self.epsilon
        )

        print("Ležaj:", self.designation)
        print("C1:", self.C1)
        print("Vrijeme:", self.Lh)

        print(f"{'' if self.C1 <= self.C else 'NE '}ZADOVOLJAVA!")
