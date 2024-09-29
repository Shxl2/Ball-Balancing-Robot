import numpy as np
import math

class StewartPlatform:
    def __init__(self, la, lb, lc, ld):
        self.la = la
        self.lb = lb
        self.lc = lc
        self.ld = ld

        self.a1 = np.array([0, la, 0])
        self.a2 = np.array([-math.sqrt(3) * la / 2, -la / 2, 0])
        self.a3 = np.array([math.sqrt(3) * la / 2, -la / 2, 0])

        self.d1 = np.array([0, ld, 0])
        self.d2 = np.array([-math.sqrt(3) * ld / 2, -ld / 2, 0])
        self.d3 = np.array([math.sqrt(3) * ld / 2, -ld / 2, 0])

    def getAngles(self, i, j, h):
        k = self.getK(i, j)

        R = self.getR(i, j, k)

        P = self.getP(i, j, k, h)
      
        f1 = np.matmul(R, self.d1) + P - self.a1
        f2 = np.matmul(R, self.d2) + P - self.a2
        f3 = np.matmul(R, self.d3) + P - self.a3

        return np.array([self.getAngle(f1, self.a1), self.getAngle(f2, self.a2), self.getAngle(f3, self.a3)])

    def getK(self, i, j):
        si, sj = math.sin(i), math.sin(j)
        ci, cj = math.cos(i), math.cos(j)

        return math.atan((si * sj) / (ci + cj))
    
    def getR(self, i, j, k):
        si, sj, sk = math.sin(i), math.sin(j), math.sin(k)
        ci, cj, ck = math.cos(i), math.cos(j), math.cos(k)

        return np.array([[cj * ck, si * sj * sk - ci * sk, ci * sj * ck + si * sk], 
                 [cj * sk, si * sj * sk + ci * ck, ci * sj * sk - si * ck], 
                 [-sj, si * cj, ci * cj]])
        
    def getP(self, i, j, k, h):
        si, sj, sk = math.sin(i), math.sin(j), math.sin(k)
        ci, cj, ck = math.cos(i), math.cos(j), math.cos(k)

        return np.array([-self.ld * cj * sk, (self.ld / 2) * (ci * ck + si * sj * sk - cj * ck), h])
    
    def getAngle(self, fi, ai):
        magfi = np.sqrt(fi.dot(fi))

        theta1 = math.acos(np.dot(ai, fi) / (self.la * magfi))
        theta2 = math.acos((self.lb * self.lb + magfi * magfi - self.lc * self.lc) / (2 * self.lb * magfi))

        return theta1 - theta2

