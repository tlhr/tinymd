#![allow(dead_code)]

extern crate linalg;

use std::thread;
use std::sync::mpsc;
use std::io::prelude::*;
use std::io::BufReader;
use std::fs::File;
use std::fs::OpenOptions;
use linalg::{Vec3, Mat3};

/// Reads a 3 x SIZE array of floats from a textfile, returning a Mat3.
fn read_input(filename: String) -> Mat3 {
    let mut f = BufReader::new(File::open(filename).unwrap());
    let mut s = String::new();
    f.read_line(&mut s).unwrap();

    let arr: Mat3 = Mat3::from_vec(f.lines().map(|l| {
        Vec3::from_vec(l.unwrap()
                        .split(char::is_whitespace)
                        .map(|s| s.trim())
                        .filter(|s| !s.is_empty())
                        .map(|number| number.parse().unwrap())
                        .collect())
    }).collect());
    return arr;
}

/// Writes data into an xyz file.
///
/// Takes a Mat3 array and appends it to an existing or new xyz file,
/// including a comment line and the atom type.
fn write_data(x: Mat3, outputfile: String, atom: String) {
    let mut f = OpenOptions::new()
                    .write(true)
                    .append(true)
                    .create(true)
                    .open(outputfile)
                    .unwrap();
    {
        write!(f, "{}\n", x.dat.len()).unwrap();
        write!(f, "# {}\n", "much data").unwrap();
        for i in 0..x.dat.len() {
            write!(f, "{} {} {} {}\n", atom,
                   x.dat[i].x, x.dat[i].y, x.dat[i].z).unwrap();
        }
    }
}

#[derive(Copy, Clone, Debug)]
struct VerletEnsemble {
    y0: Mat3,
    z0: Mat3,
    y: Mat3,
    z: Mat3,
    dt: f64,
}

impl VerletEnsemble {
    fn new(y0: Mat3, z0: Mat3, y: Mat3, z: Mat3, dt: f64) -> VerletEnsemble {
        VerletEnsemble {y0: y0, z0: z0, y: y, z: z, dt: dt}
    }

    fn init(self) -> Vec<thread::JoinHandle<_>> {
        (0..self.y0.dat.len()).map(|atom| {
            let (mut fyv, mut yatom, mut zatom): (Vec3, Vec3, Vec3);
            thread::spawn(move || {
                fyv = self.calc_fy(atom);
                yatom = self.calc_y(fyv, atom);
                zatom = self.calc_z(fyv, atom);
            })
        }).collect()
    }

    fn calc_fy(self, i: usize) -> Vec3 {
        lj(&self.y0, &i)
    }

    fn calc_y(self, fyv: Vec3, i: usize) -> Vec3 {
        self.z[i] * self.dt + self.y0[i] + fyv * self.dt * 0.5
    }

    fn calc_z(self, fyv: Vec3, i: usize) -> Vec3 {
        self.z0[i] + fyv * self.dt * 0.5 + lj(&self.y, &i) * self.dt * 0.5
    }
}

/// Verlet integrator for two coupled first-order diff. equations.
fn verlet(y0: Mat3, z0: Mat3, dt: f64, f: fn(&Mat3, &usize) -> Vec3) -> (Mat3, Mat3) {
    let (mut y, mut z, mut fym) = (Mat3::zeros(), Mat3::zeros(), Mat3::zeros());
    let mut fyv: Vec3;

    for i in 0..y0.dat.len() {
        fyv = f(&y0, &i);
        y[i] = y0[i] + z[i] * dt + fyv * dt * 0.5;
        fym[i] = fyv;
    }

    for i in 0..y0.dat.len() {
        z[i] = z0[i] + fym[i] * dt * 0.5 + f(&y, &i) * dt * 0.5;
    }
    return (y, z);
}

/// Runge-Kutta 4th order integrator for two coupled first-order
/// differential equations.
///
/// `fn1`: momentum function
/// `fn2`: potential function
fn rk4(y0: Mat3,
       z0: Mat3,
       dt: f64,
       fn1: fn(Mat3, usize) -> Vec3,
       fn2: fn(Mat3, usize) -> Vec3)
       -> (Mat3, Mat3) {

    let (mut y, mut z): (Mat3, Mat3) = (Mat3::zeros(), Mat3::zeros());
    let (mut ktemp, mut ltemp): (Vec3, Vec3);
    let (mut k1, mut k2, mut k3, mut k4): (Vec3, Vec3, Vec3, Vec3);
    let (mut l1, mut l2, mut l3, mut l4): (Vec3, Vec3, Vec3, Vec3);
    let (dt2, dt3, dt6): (f64, f64, f64) = (dt / 2., dt / 3., dt / 6.);

    for i in 0..y0.dat.len() {
        ktemp = y0[i];
        ltemp = z0[i];
        k1 = fn1(z0, i);
        l1 = fn2(y0, i);
        k2 = fn1(z0 + l1 * dt2, i);
        l2 = fn2(y0 + k1 * dt2, i);
        k3 = fn1(z0 + l2 * dt2, i);
        l3 = fn2(y0 + k2 * dt2, i);
        k4 = fn1(z0 + l3 * dt, i);
        l4 = fn2(y0 + k3 * dt, i);
        ktemp = ktemp + k1 * dt6 + k2 * dt3 + k3 * dt3 + k4 * dt6;
        ltemp = ltemp + l1 * dt6 + l2 * dt3 + l3 * dt3 + l4 * dt6;
        y.dat[i] = ktemp;
        z.dat[i] = ltemp;
    }
    return (y, z);
}

/// Calculates the atom energy of the current state
/// given the coordinate and impulse matrix and index i.
fn energy(i: usize, q: Mat3, p: Mat3) -> f64 {
    let mut v: f64 = 0.;
    let mut r: f64;
    let sigma: f64 = 1.;
    let e: f64 = 1.;
    let m: f64 = 10e5;

    for j in 0..q.dat.len() {
        if i != j {
            r = (q.dat[i] - q.dat[j]).norm();
            v = v + (-4. * e * ((sigma / r).powi(12) - (sigma / r).powi(6)));
        }
    }
    return v + p[i].dot(p[i]) / (2. * m);
}

/// Calculates the updated momentum at position i.
fn mom(p: Mat3, i: usize) -> Vec3 {
    let m: f64 = 10e5;
    return p[i] / m;
}

/// Lennard-Jones potential function for rare-gas clusters.
fn lj(q: &Mat3, j: &usize) -> Vec3 {
    let m: f64 = 10e5;
    let mut pot: Vec3 = Vec3::new(0., 0., 0.);
    let mut r: f64;
    let mut diff: Vec3;

    for i in 0..q.dat.len() {
        if &i != j {
            diff = q[j] - q[i];
            r = diff.norm();
            pot = pot + (diff * -48. / r.powi(14) + diff * 24. / r.powi(8));
        }
    }
    return -pot / m;
}

/// Converts an impulse matrix to a velocity matrix.
fn p_to_v(p: Mat3) -> Mat3 {
    let m: f64 = 10e5;
    let mass: f64 = m / p.dat.len() as f64;
    return p / mass;
}

/// Converts a velocity matrix to an impulse matrix.
fn v_to_p(v: Mat3) -> Mat3 {
    let m: f64 = 10e5;
    let mass: f64 = m / v.dat.len() as f64;
    return v * mass;
}

/// Calculates the root-mean-square velocity of the cluster.
fn v_rms(v: Mat3) -> f64 {
    let mut vlist = vec![];
    let mut vtot: f64 = 0.;
    for i in 0..v.dat.len() {
        vlist.push(v[i].norm());
    }
    for i in 0..vlist.len() {
        vtot = vtot + vlist[i].powi(2);
    }
    return (vtot / v.dat.len() as f64).sqrt();
}

/// Calculates the center-of-mass velocity of the cluster.
fn v_cm(v: Mat3) -> Vec3 {
    let mut vec: Vec3 = Vec3::new(0., 0., 0.);
    for i in 0..v.dat.len() {
        vec = vec + v[i];
    }
    return vec / v.dat.len() as f64;
}

/// Calculates the total energy of the cluster.
fn etot(x: Mat3, p: Mat3) -> f64 {
    let mut e: f64 = 0.;
    for i in 0..p.dat.len() {
        e = e + energy(i, x, p);
    }
    return e;
}

/// TODO Incorrect results!
fn pmd(xinit: Mat3, pinit: Mat3, outputfile: &str, atom: &str,
       steps: u64, dt: f64, res: u64, verbose: bool) {
    let (mut x, mut v) = (xinit, p_to_v(pinit));
    let (mut x0, mut v0) = (x, v);
    let mut handles: Vec<thread::JoinHandle<_>> = vec![];
    let size = x.dat.len();

    let (tx, rx) = mpsc::channel();
    let (tv, rv) = mpsc::channel();

    for atom in 0..size {
        let (tx, tv) = (tx.clone(), tv.clone());
        handles.push(thread::spawn(move || {
            let mut fxv: Vec3;
            let (mut xatom, mut vatom): (Vec3, Vec3);
            for step in 0..steps {
                fxv = lj(&x0, &atom);
                xatom = x0[atom] + v[atom] * dt + fxv * dt * 0.5;
                tx.send(xatom).unwrap();
                thread::park();
                vatom = v0[atom] + fxv * dt * 0.5 + lj(&x, &atom) * dt * 0.5;
                tv.send(vatom).unwrap();
                thread::park();
            }
        }));
    }

    let master = thread::spawn(move || {
        for step in 0..steps {
            let (mut xtemp, mut vtemp) = (vec![], vec![]);
            for _ in 0..size {
                xtemp.push(rx.recv().unwrap());
            }
            x = Mat3::from_vec(xtemp);
            x0 = x;

            for h in &handles {
                h.thread().unpark();
            }

            for _ in 0..size {
                vtemp.push(rv.recv().unwrap());
            }
            v = Mat3::from_vec(vtemp);
            v = v - v_cm(v);
            v0 = v;

            for h in &handles {
                h.thread().unpark();
            }
            if verbose {
                println!("Round: {}/{} :: v_rms: {} :: E: {}",
                         step + 1, steps, v_rms(v), etot(x, v_to_p(v)));
            }
            // write data with desired resolution as xyz file
            // if step % res == 0 {
            //     write_data(x, outputfile.to_owned(), atom.to_owned());
            // }
        }
        // Join all threads to end simulation
        for h in handles {
            h.join().unwrap();
        }
    });
    master.join().unwrap();
}

/// Runs the molecular dynamics simulation
fn md(xinit: Mat3, pinit: Mat3, outputfile: &str, atom: &str,
      steps: u64, dt: f64, res: u64, verbose: bool) {
    let (mut x, mut v) = (xinit, p_to_v(pinit));
    for t in 0..steps {
        let (xi, vi) = verlet(x, v, dt, lj);
        x = xi;
        v = vi;
        // center-of-mass velocity correction
        v = v - v_cm(v);
        if verbose {
            println!("Round: {}/{} :: v_rms: {} :: E: {}",
                     t + 1, steps, v_rms(v), etot(x, v_to_p(v)));
        }
        // write data with desired resolution as xyz file
        if t % res == 0 {
            write_data(x, outputfile.to_owned(), atom.to_owned());
        }
    }
}

fn main() {
    let x = read_input("input/lj20.dat".to_owned());
    let mut v = Mat3::zeros();
    v[1] = Vec3::new(0.00, 0.00, 0.05);
    let p: Mat3 = v_to_p(v);
    md(x, p, "out/test.xyz", "Ar", 10000, 1., 10, true);
}
