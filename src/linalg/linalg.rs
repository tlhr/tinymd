use std::ops::{Add, Sub, Mul, Div, Neg, Index, IndexMut};

pub const SIZE: usize = 20;

#[derive(Copy, Clone, Debug)]
pub struct Vec3 {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Vec3 {
    pub fn new(x: f64, y: f64, z: f64) -> Vec3 {
        Vec3 { x: x, y: y, z: z }
    }

    pub fn from_vec(vec: Vec<f64>) -> Vec3 {
        Vec3 { x: vec[0], y: vec[1], z: vec[2] }
    }

    pub fn zeros() -> Vec3 {
        Vec3 { x: 0., y: 0., z: 0. }
    }

    pub fn dot(&self, other: Vec3) -> f64 {
        return self.x * other.x +
               self.y * other.y +
               self.z * other.z
    }

    pub fn norm(&self) -> f64 {
        let norm: f64 = self.x * self.x +
                        self.y * self.y +
                        self.z * self.z;
        return norm.sqrt()
    }
}

impl Neg for Vec3 {
    type Output = Vec3;

    fn neg(self) -> Vec3 {
        Vec3 { x: self.x * -1.,
               y: self.y * -1.,
               z: self.z * -1. }
    }
}

impl Add for Vec3 {
    type Output = Vec3;

    fn add(self, other: Vec3) -> Vec3 {
        Vec3 { x: self.x + other.x,
               y: self.y + other.y,
               z: self.z + other.z }
    }
}

impl Add<f64> for Vec3 {
    type Output = Vec3;

    fn add(self, other: f64) -> Vec3 {
        Vec3 { x: self.x + other,
               y: self.y + other,
               z: self.z + other }
    }
}

impl<'b> Add for &'b Vec3 {
    type Output = Vec3;

    fn add(self, other: &'b Vec3) -> Vec3 {
        Vec3 { x: self.x + other.x,
               y: self.y + other.y,
               z: self.z + other.z }
    }
}

impl<'a> Add<&'a Vec3> for Vec3 {
    type Output = Vec3;

    fn add(self, other: &'a Vec3) -> Vec3 {
        Vec3 { x: self.x + other.x,
               y: self.y + other.y,
               z: self.z + other.z }
    }
}

impl Sub for Vec3 {
    type Output = Vec3;

    fn sub(self, other: Vec3) -> Vec3 {
        Vec3 { x: self.x - other.x,
               y: self.y - other.y,
               z: self.z - other.z }
    }
}

impl Sub<f64> for Vec3 {
    type Output = Vec3;

    fn sub(self, other: f64) -> Vec3 {
        Vec3 { x: self.x - other,
               y: self.y - other,
               z: self.z - other }
    }
}

impl Mul for Vec3 {
    type Output = Vec3;

    fn mul(self, other: Vec3) -> Vec3 {
        Vec3 { x: self.x * other.x,
               y: self.y * other.y,
               z: self.z * other.z }
    }
}

impl Mul<f64> for Vec3 {
    type Output = Vec3;

    fn mul(self, other: f64) -> Vec3 {
        Vec3 { x: self.x * other,
               y: self.y * other,
               z: self.z * other }
    }
}

impl<'b> Mul<f64> for &'b Vec3 {
    type Output = Vec3;

    fn mul(self, other: f64) -> Vec3 {
        Vec3 { x: self.x * other,
               y: self.y * other,
               z: self.z * other }
    }
}

impl<'b> Mul<&'b f64> for Vec3 {
    type Output = Vec3;

    fn mul(self, other: &f64) -> Vec3 {
        Vec3 { x: self.x * other,
               y: self.y * other,
               z: self.z * other }
    }
}

impl<'b> Mul<&'b f64> for &'b Vec3 {
    type Output = Vec3;

    fn mul(self, other: &f64) -> Vec3 {
        Vec3 { x: self.x * other,
               y: self.y * other,
               z: self.z * other }
    }
}

impl<'b> Mul for &'b Vec3 {
    type Output = Vec3;

    fn mul(self, other: &'b Vec3) -> Vec3 {
        Vec3 { x: self.x * other.x,
               y: self.y * other.y,
               z: self.z * other.z }
    }
}

impl Div for Vec3 {
    type Output = Vec3;

    fn div(self, other: Vec3) -> Vec3 {
        Vec3 { x: self.x / other.x,
               y: self.y / other.y,
               z: self.z / other.z }
    }
}

impl Div<f64> for Vec3 {
    type Output = Vec3;

    fn div(self, other: f64) -> Vec3 {
        Vec3 { x: self.x / other,
               y: self.y / other,
               z: self.z / other }
    }
}

#[derive(Copy, Clone, Debug)]
pub struct Mat3 {
    pub dat: [Vec3; SIZE],
}

impl Mat3 {
    pub fn new(veclist: [Vec3; SIZE]) -> Mat3 {
        Mat3 { dat: veclist }
    }

    pub fn from_vec(vec: Vec<Vec3>) -> Mat3 {
        let mut mat = Mat3::zeros();
        for i in 0..vec.len() {
            mat.dat[i] = vec[i];
        }
        return mat
    }

    pub fn zeros() -> Mat3 {
        return Mat3::new([Vec3::new(0., 0., 0.); SIZE]);
    }
}

impl Add<Vec3> for Mat3 {
    type Output = Mat3;

    fn add(self, other: Vec3) -> Mat3 {
        let mut mat = Mat3::zeros();
        for i in 0..self.dat.len() {
            mat.dat[i] = self.dat[i] + other;
        }
        return mat
    }
}

impl Sub<Vec3> for Mat3 {
    type Output = Mat3;

    fn sub(self, other: Vec3) -> Mat3 {
        let mut mat = Mat3::zeros();
        for i in 0..self.dat.len() {
            mat.dat[i] = self.dat[i] - other;
        }
        return mat
    }
}

impl Mul<f64> for Mat3 {
    type Output = Mat3;

    fn mul(self, other: f64) -> Mat3 {
        let mut mat = Mat3::zeros();
        for i in 0..self.dat.len() {
            mat.dat[i] = self.dat[i] * other;
        }
        return mat
    }
}

impl Div<f64> for Mat3 {
    type Output = Mat3;

    fn div(self, other: f64) -> Mat3 {
        let mut mat = Mat3::zeros();
        for i in 0..self.dat.len() {
            mat.dat[i] = self.dat[i] / other;
        }
        return mat
    }
}

impl Index<usize> for Mat3 {
    type Output = Vec3;

    fn index<'a>(&'a self, _index: usize) -> &'a Vec3 {
        let opt = self.dat.get(_index);
        return opt.expect("Matrix index out of bounds!")
    }
}

impl<'b> Index<&'b usize> for Mat3 {
    type Output = Vec3;

    fn index<'a>(&'a self, _index: &usize) -> &'a Vec3 {
        let i = *_index;
        let opt = self.dat.get(i);
        return opt.expect("Matrix index out of bounds!")
    }
}

impl IndexMut<usize> for Mat3 {
    fn index_mut<'a>(&'a mut self, _index: usize) -> &'a mut Vec3 {
        let opt = self.dat.get_mut(_index);
        return opt.expect("Matrix index out of bounds!")
    }
}
