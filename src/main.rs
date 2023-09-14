use plotters::prelude::*;
use rand::prelude::*;
use rand_xoshiro::rand_core::SeedableRng;
use rand_xoshiro::Xoshiro256Plus;
use std::time::Instant;

/*
   This code is largely based off of an implementation I fist saw written in
   julia's Ising2D.jl written by genkuroki: https://github.com/genkuroki/Ising2D.jl/blob/master/src/Ising2D.jl
   I wanted to see how a fairly naive port to rust would do as most the work
   is done in "low level" ways using arrays and for loops. The performance of
   this rust code is lackluster, taking twice the time of the julia implementation.
   My rewritten stripped down julia implementation is included along side this one
   for anyone iterested in testing for the comparison.
   */

fn main() {
    let t = 2. / f64::ln(1. + 2f64.sqrt()); //The critical temperature of a 2D Ising model.

    run(true, t);
}
/*
   List of constants that will apply for the rest of the program
   J is the coupling constant of the interation, it can be a float
   in principle, but I will use 1 for the classic ferromagnetic case.
   */
const J: i8 = 1;
const STEPS: usize = 1000;
const SIDE: usize = 1000; // Making a default for square arrays
const NPIXELS: u32 = SIDE as u32; // Used for giving the size of a side of the PNG.
const NROWS: usize = SIDE;
const NCOLUMNS: usize = SIDE;
const LEN: usize = NROWS * NCOLUMNS;

/*
   Defining a function to run the simulation inside of main. Order
   determines if it is a "hot" or "cold" initial state.
   */
fn run(order: bool, t: f64) {
    /*
       initialize the array we will be using through the rest of the program.
       The rng is needed for determining if the state flips or not and for
       randomizing the initial state when order is passed as true.
       */
    let mut arr = [0i8; LEN];
    let mut rng = Xoshiro256Plus::from_entropy();

    let beta = 1. / t; // beta is a convenience variable for inverse Temp.

    /*
       When order is true the array is randomized between up (1) and down (-1)
       spin states. Otherwise the spins all point up.
       */

    if order {
        for site in arr.iter_mut() {
            *site = if rng.gen_bool(0.5) { 1i8 } else { -1i8 }
        }
    } else {
        arr = [1i8; LEN]
    }

    // Show the image before the interations for comparison.
    plot(&arr, String::from("before.png").as_str()).unwrap();

    let start = Instant::now();
    /*
       Create a static array of probabilites based on the local energy of a given
       site. This can be done because the local energy is one of a discrete number
       of possible outcomes and enables one to check the energy of the site as
       opposed to the whole lattice when determining whether or not to flip the spin.
       */
    // We move it to u64 on suggestion as comparing u64s is faster than floats.
    let mut probs = [0u64; 9];
    let mut increment: f64 = -4.0;
    for prob in &mut probs {
        let ptemp = f64::exp(-2. * beta * increment);
        *prob = (2f64.powi(64) * ptemp) as u64; //This is essentially witchcraft.
        increment += 1.;
    }// The generating of the 9 random numbers should be possible to make parallel
    // TODO: Try generating several random floats in parallel.

    /*
       This is the major loop of the model. STEPS is the number of times we iterate through
       the lattice trying to flip each spin. j is the iterating over the "columns"
       of the array, and i is iterating over the "rows", though the array is technically
       1D. First the energy of each site is calculated based on the Ising model's
       hamiltonian, h = site * (nn + ss + ee + ww), where h is the energy of the site,
       site is the spin of the value of the site in question, and nn, ss, ee, ww are
       the site above, below, right and left of the current position respectively.
       I'm using wrapping boundary conditions to mitigate edge effects, meaning I
       am treating the right edge of the 2D array as touching the left edge, and
       the top wraps to the bottom. If the energy is lowered by fliping the site's spin we
       do so, otherwise we flip the spin with probability exp(-2 * energy at the site * beta).
       */


    for _ in 0..STEPS {
        for i in 0..NROWS {
            for j in 0..NCOLUMNS {
                // precalculate i's for the current iteration
                let inorth = ((i + 1) % NROWS) * NROWS;
                let isouth = if i == 0 {
                    (NROWS - 1) * NROWS
                } else {
                    (i - 1) * NROWS
                };
                let i = i * NROWS; //shadow i with its current row value.

                // assert! helps elide bounds checks
                assert!(i + j < arr.len());

                let jeast = (j + 1) % NCOLUMNS;
                let jwest = if j == 0 {
                    (j + NCOLUMNS - 1) % NCOLUMNS
                } else {
                    j - 1
                };

                let nn = &arr[inorth + j];
                let ss = &arr[isouth + j];
                let ee = &arr[i + jeast];
                let ww = &arr[i + jwest];
                let site = &arr[i + j];

                let en = J * site * (nn + ss + ww + ee);
                let pcomp = &rng.gen::<u64>();

                let k = 4 + en;
                assert!((k as usize) < probs.len());
                let flip = *pcomp < probs[k as usize];

                arr[i + j] = if flip { -arr[i + j] } else { arr[i + j] };
            }
        }
    }

    // Check how long the program took to run.
    let elapsed = Instant::now() - start;
    println!("the whole program took {:#?} seconds to run.", elapsed);

    // Plot the final state the system is in.
    plot(&arr, String::from("after.png").as_str()).unwrap();
}

/*
   Defining a basic plotting function for the array as a png. This maps the array to a 2D
   histogram where an up spin (1) is white, and a down spin (-1) is teal.
   */
fn plot(&arr: &[i8; LEN], name: &str) -> Result<(), Box<dyn std::error::Error>> {
    let root_drawing_area = BitMapBackend::new(name, (NPIXELS, NPIXELS)).into_drawing_area();

    let child_drawing_areas = root_drawing_area.split_evenly((NROWS, NCOLUMNS));

    for (area, i) in child_drawing_areas.into_iter().zip(0..LEN) {
        if arr[i] == 1i8 {
            area.fill(&WHITE)?;
        } else {
            let teal = RGBColor(0u8, 128u8, 128u8);
            area.fill(&teal)?;
        }
    }

    Ok(
        ()
      )
}
