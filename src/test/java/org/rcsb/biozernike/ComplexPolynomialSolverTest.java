package org.rcsb.biozernike;

import org.junit.jupiter.api.Test;
import org.rcsb.biozernike.complex.Complex;
import org.rcsb.biozernike.complex.PolynomialSolver;

import java.util.Arrays;
import java.util.List;

import static org.junit.jupiter.api.Assertions.*;

public class ComplexPolynomialSolverTest {

	@Test
	public void roots() {
		double eps = 1E-10;
		PolynomialSolver solver = new PolynomialSolver();

		// simple quadratic
		Complex[] coeffs = {
				new Complex(1,0),
				new Complex(1,0),
				new Complex(1,0)
		};

		Complex[] roots = solver.roots(coeffs);
		assertEquals(2, roots.length);
		assertEquals(-0.5, roots[0].getReal(), eps);
		assertEquals(-0.5, roots[1].getReal(), eps);
		assertEquals(Math.sqrt(3)/2.0, Math.abs(roots[0].getImaginary()), eps);
		assertEquals(roots[0].getImaginary(), -roots[1].getImaginary(), eps);

		// roots have imaginary parts only
		coeffs = new Complex[]{
				new Complex(1,0),
				new Complex(0,-5),
				new Complex(-6,0)
		};

		roots = solver.roots(coeffs);
		assertEquals(2, roots.length);
		assertEquals(0,roots[0].getReal(),eps);
		assertEquals(0,roots[1].getReal(),eps);
		assertEquals(3,roots[0].getImaginary(),eps);
		assertEquals(2,roots[1].getImaginary(),eps);

		// roots have real parts only
		coeffs = new Complex[]{
				new Complex(1,0),
				new Complex(5,0),
				new Complex(0,0),
				new Complex(-20,0),
				new Complex(-10,0),
				new Complex(2,0)
		};

		List<Double> expectedRootsRe = Arrays.asList(
				1.8937458708571862,
				-2.5125791239422717,
				-0.7185948800821463,
				0.15328930841789012,
				-3.8158611752506584
		);

		roots = solver.roots(coeffs);
		assertEquals(5, roots.length);
		for (Complex root: roots) {
			assertEquals(0,root.getImaginary(),eps);
			assertTrue(expectedRootsRe.contains(root.getReal()));
		}

		// mixed
		coeffs = new Complex[]{
				new Complex(1,1),
				new Complex(2,-1),
				new Complex(3,1),
				new Complex(4,-1),
				new Complex(5,1)
		};

		List<Complex> expectedRoots = Arrays.asList(
				new Complex(0.56637821781158,1.8328367580651),
				new Complex(0.28462734463903,-1.2136721483554),
				new Complex(-0.60124476269297,1.4663824544435),
				new Complex(-0.74976079975764,-0.58554706415324)
		);

		roots = solver.roots(coeffs);
		assertEquals(4, roots.length);
		for (Complex root: roots) {
			boolean rootFound = false;
			for (Complex expectedRoot:expectedRoots) {
				if (Math.abs(root.getReal()-expectedRoot.getReal())<eps && Math.abs(root.getImaginary()-expectedRoot.getImaginary())<eps) {
					rootFound = true;
					break;
				}
			}
			assertTrue(rootFound);
		}

	}
}