extern crate mpolynomial;

#[cfg(test)]
mod tests {
    #[test]
    fn binomail_2() {
        // Test  if (x+y)^2 = x^2 + 2 x y + y ^2
        let n_var = 2;
        let n = 2;
        let mut mpoly_exact = mpolynomial::MPolynomial::<f64>::new(n_var);
        mpoly_exact.add(&vec![2, 0], 1.0);
        mpoly_exact.add(&vec![1, 1], 2.0);
        mpoly_exact.add(&vec![0, 2], 1.0);
        // using power function
        let mpoly_test = mpolynomial::MPolynomial::linear_pown(&[1.0, 1.0], &[1, 2], n_var, n);
        assert_eq!(mpoly_exact, mpoly_test);
    }
    #[test]
    fn binomail_3() {
        // Test  if (x+y)^3 = x^3 + 3 x^2 y + 3 x y^2 + y ^2
        let n_var = 2;
        let n = 3;
        let mut mpoly_exact = mpolynomial::MPolynomial::<f64>::new(n_var);
        mpoly_exact.add(&vec![3, 0], 1.0);
        mpoly_exact.add(&vec![2, 1], 3.0);
        mpoly_exact.add(&vec![1, 2], 3.0);
        mpoly_exact.add(&vec![0, 3], 1.0);
        // using power function
        let mpoly_test = mpolynomial::MPolynomial::linear_pown(&[1.0, 1.0], &[1, 2], n_var, n);
        assert_eq!(mpoly_exact, mpoly_test);
    }
    #[test]
    fn multinomial_3_2() {
        // Test  if (x1+x2+x3)^2
        let n_var = 3;
        let n = 2;
        let mut mpoly_exact = mpolynomial::MPolynomial::<f64>::new(n_var);
        mpoly_exact.add(&vec![2, 0, 0], 1.0);
        mpoly_exact.add(&vec![0, 2, 0], 1.0);
        mpoly_exact.add(&vec![0, 0, 2], 1.0);
        mpoly_exact.add(&vec![1, 1, 0], 2.0);
        mpoly_exact.add(&vec![0, 1, 1], 2.0);
        mpoly_exact.add(&vec![1, 0, 1], 2.0);
        // using power function
        let mpoly_test =
            mpolynomial::MPolynomial::linear_pown(&[1.0, 1.0, 1.0], &[1, 2, 3], n_var, n);
        assert_eq!(mpoly_exact, mpoly_test);
    }
    #[test]
    fn multinomial_3_3() {
        // Test  if (x1+x2+x3)^3
        let n_var = 3;
        let n = 3;
        let mut mpoly_exact = mpolynomial::MPolynomial::<f64>::new(n_var);
        mpoly_exact.add(&vec![3, 0, 0], 1.0);
        mpoly_exact.add(&vec![2, 1, 0], 3.0);
        mpoly_exact.add(&vec![2, 0, 1], 3.0);
        mpoly_exact.add(&vec![1, 2, 0], 3.0);
        mpoly_exact.add(&vec![1, 0, 2], 3.0);
        mpoly_exact.add(&vec![1, 1, 1], 6.0);
        mpoly_exact.add(&vec![0, 3, 0], 1.0);
        mpoly_exact.add(&vec![0, 0, 3], 1.0);
        mpoly_exact.add(&vec![0, 2, 1], 3.0);
        mpoly_exact.add(&vec![0, 1, 2], 3.0);
        // using power function
        let mpoly_test =
            mpolynomial::MPolynomial::linear_pown(&[1.0, 1.0, 1.0], &[3, 1, 2], n_var, n);
        assert_eq!(mpoly_exact, mpoly_test);
    }

    #[test]
    fn replace_3() {
        // test rotation on (1+3x1+5x2+7x3)^6 with x1 -> 11 - 13x1 + 17x2 + 19x3
        let n_var = 3;
        let n = 6;
        let mut mpoly_exact = mpolynomial::MPolynomial::<f64>::new(n_var);
        mpoly_exact.add(&vec![0, 0, 0], 1544804416.0);
        mpoly_exact.add(&vec![1, 0, 0], -10631889216.0); // x1
        mpoly_exact.add(&vec![2, 0, 0], 30488505840.0); // x1^2
        mpoly_exact.add(&vec![3, 0, 0], -46629479520.0); //x1^3
        mpoly_exact.add(&vec![4, 0, 0], 40115066940.0); //x1^4
        mpoly_exact.add(&vec![5, 0, 0], -18405736596.0); //x1^5
        mpoly_exact.add(&vec![6, 0, 0], 3518743761.0); //x1^6
        mpoly_exact.add(&vec![0, 1, 0], 15266302464.0); //x2
        mpoly_exact.add(&vec![1, 1, 0], -87556734720.0); //x1 x2
        mpoly_exact.add(&vec![2, 1, 0], 200865450240.0); //x1^2 x2
        mpoly_exact.add(&vec![3, 1, 0], -230404487040.0); //x1^3 x2
        mpoly_exact.add(&vec![4, 1, 0], 132143749920.0); //x1^4 x2
        mpoly_exact.add(&vec![5, 1, 0], -30315330864.0); //x1^5 x2
        mpoly_exact.add(&vec![0, 2, 0], 62861245440.0); //x2^2
        mpoly_exact.add(&vec![1, 2, 0], -288422184960.0); //x1 x2^2
        mpoly_exact.add(&vec![2, 2, 0], 496255818240.0); //x1^2 x2^2
        mpoly_exact.add(&vec![3, 2, 0], -379489743360.0); //x1^3 x2^2
        mpoly_exact.add(&vec![4, 2, 0], 108824264640.0); //x1^4 x2^2
        mpoly_exact.add(&vec![0, 3, 0], 138048225280.0); //x2^3
        mpoly_exact.add(&vec![1, 3, 0], -475048304640.0); //x1 x2^3
        mpoly_exact.add(&vec![2, 3, 0], 544908349440.0); //x1^2 x2^3
        mpoly_exact.add(&vec![3, 3, 0], -208347310080.0); //x1^3 x2^3
        mpoly_exact.add(&vec![0, 4, 0], 170530160640.0); //x2^4
        mpoly_exact.add(&vec![1, 4, 0], -391216250880.0); //x1 x2^4
        mpoly_exact.add(&vec![2, 4, 0], 224374026240.0); //x1^2 x2^4
        mpoly_exact.add(&vec![0, 5, 0], 112349282304.0); //x2^5
        mpoly_exact.add(&vec![1, 5, 0], -128871235584.0); //x1 x2^5
        mpoly_exact.add(&vec![0, 6, 0], 30840979456.0); //x2^6
        mpoly_exact.add(&vec![0, 0, 1], 17447202816.0); //x3
        mpoly_exact.add(&vec![1, 0, 1], -100064839680.0); //x1 x3
        mpoly_exact.add(&vec![2, 0, 1], 229560514560.0); //x1^2 x3
        mpoly_exact.add(&vec![3, 0, 1], -263319413760.0); //x1^3 x3
        mpoly_exact.add(&vec![4, 0, 1], 151021428480.0); //x1^4 x3
        mpoly_exact.add(&vec![5, 0, 1], -34646092416.0); //x1^5 x3
        mpoly_exact.add(&vec![0, 1, 1], 143682846720.0); //x2 x3
        mpoly_exact.add(&vec![1, 1, 1], -659250708480.0); //x1 x2 x3
        mpoly_exact.add(&vec![2, 1, 1], 1134299013120.0); //x1^2 x2 x3
        mpoly_exact.add(&vec![3, 1, 1], -867405127680.0); //x1^3 x2 x3
        mpoly_exact.add(&vec![4, 1, 1], 248741176320.0); //x1^4 x2 x3
        mpoly_exact.add(&vec![0, 2, 1], 473308200960.0); //x2^2 x3
        mpoly_exact.add(&vec![1, 2, 1], -1628737044480.0); //x1 x2^2 x3
        mpoly_exact.add(&vec![2, 2, 1], 1868257198080.0); //x1^2 x2^2 x3
        mpoly_exact.add(&vec![3, 2, 1], -714333634560.0); //x1^3 x2^2 x3
        mpoly_exact.add(&vec![0, 3, 1], 779566448640.0); //x2^3 x3
        mpoly_exact.add(&vec![1, 3, 1], -1788417146880.0); //x1 x2^3 x3
        mpoly_exact.add(&vec![2, 3, 1], 1025709834240.0); //x1^2 x2^3 x3
        mpoly_exact.add(&vec![0, 4, 1], 641995898880.0); //x2^4 x3
        mpoly_exact.add(&vec![1, 4, 1], -736407060480.0); //x1 x2^4 x3
        mpoly_exact.add(&vec![0, 5, 1], 211481001984.0); //x2^5 x3
        mpoly_exact.add(&vec![0, 0, 2], 82104483840.0); //x3^2
        mpoly_exact.add(&vec![1, 0, 2], -376714690560.0); //x1 x3^2
        mpoly_exact.add(&vec![2, 0, 2], 648170864640.0); //x1^2 x3^2
        mpoly_exact.add(&vec![3, 0, 2], -495660072960.0); //x1^3 x3^2
        mpoly_exact.add(&vec![4, 0, 2], 142137815040.0); //x1^4 x3^2
        mpoly_exact.add(&vec![0, 1, 2], 540923658240.0); //x2 x3^2
        mpoly_exact.add(&vec![1, 1, 2], -1861413765120.0); //x1 x2 x3^2
        mpoly_exact.add(&vec![2, 1, 2], 2135151083520.0); //x1^2 x2 x3^2
        mpoly_exact.add(&vec![3, 1, 2], -816381296640.0); //x1^3 x2 x3^2
        mpoly_exact.add(&vec![0, 2, 2], 1336399626240.0); //x2^2 x3^2
        mpoly_exact.add(&vec![1, 2, 2], -3065857966080.0); //x1 x2^2 x3^2
        mpoly_exact.add(&vec![2, 2, 2], 1758359715840.0); //x1^2 x2^2 x3^2
        mpoly_exact.add(&vec![0, 3, 2], 1467419197440.0); //x2^3 x3^2
        mpoly_exact.add(&vec![1, 3, 2], -1683216138240.0); //x1 x2^3 x3^2
        mpoly_exact.add(&vec![0, 4, 2], 604231434240.0); //x2^4 x3^2
        mpoly_exact.add(&vec![0, 0, 3], 206066155520.0); //x3^3
        mpoly_exact.add(&vec![1, 0, 3], -709110005760.0); //x1 x3^3
        mpoly_exact.add(&vec![2, 0, 3], 813390888960.0); //x1^2 x3^3
        mpoly_exact.add(&vec![3, 0, 3], -311002398720.0); //x1^3 x3^3
        mpoly_exact.add(&vec![0, 1, 3], 1018209239040.0); //x2 x3^3
        mpoly_exact.add(&vec![1, 1, 3], -2335891783680.0); //x1 x2 x3^3
        mpoly_exact.add(&vec![2, 1, 3], 1339702640640.0); //x1^2 x2 x3^3
        mpoly_exact.add(&vec![0, 2, 3], 1677050511360.0); //x2^2 x3^3
        mpoly_exact.add(&vec![1, 2, 3], -1923675586560.0); //x1 x2^2 x3^3
        mpoly_exact.add(&vec![0, 3, 3], 920733614080.0); //x2^3 x3^3
        mpoly_exact.add(&vec![0, 0, 4], 290916925440.0); //x3^4
        mpoly_exact.add(&vec![1, 0, 4], -667397652480.0); //x1 x3^4
        mpoly_exact.add(&vec![2, 0, 4], 382772183040.0); //x1^2 x3^4
        mpoly_exact.add(&vec![0, 1, 4], 958314577920.0); //x2 x3^4
        mpoly_exact.add(&vec![1, 1, 4], -1099243192320.0); //x1 x2 x3^4
        mpoly_exact.add(&vec![0, 2, 4], 789200240640.0); //x2^2 x3^4
        mpoly_exact.add(&vec![0, 0, 5], 219043332096.0); //x3^5
        mpoly_exact.add(&vec![1, 0, 5], -251255586816.0); //x1 x3^5
        mpoly_exact.add(&vec![0, 1, 5], 360777252864.0); //x2 x3^5
        mpoly_exact.add(&vec![0, 0, 6], 68719476736.0); //x3^6

        // first create the expression
        let mut mpoly_test =
            mpolynomial::MPolynomial::linear_pown(&[1.0, 3.0, 5.0, 7.0], &[0, 1, 2, 3], n_var, n);
        // Now replace
        mpoly_test.replace(1, &[11.0, -13.0, 17.0, 19.0], &[0, 1, 2, 3]);
        assert_eq!(mpoly_exact, mpoly_test);
    }
    #[test]
    fn multiplication_as_power() {
        // compare (7-3x1+5x2)^6 using MPolynomial::mult and MPolynomial::linear_pown
        let n_var = 2;
        let n = 6;
        // Create two copies of (3x1 +5x2)
        let mut mpoly1 = mpolynomial::MPolynomial::new(n_var);
        mpoly1.add(&vec![0, 0], 7.0);
        mpoly1.add(&vec![1, 0], -3.0);
        mpoly1.add(&vec![0, 1], 5.0);
        let mut mpoly2 = mpolynomial::MPolynomial::new(n_var);
        mpoly2.add(&vec![0, 0], 7.0);
        mpoly2.add(&vec![1, 0], -3.0);
        mpoly2.add(&vec![0, 1], 5.0);

        // multiply
        for _ in 1..n {
            mpoly1.mult(&mpoly2);
        }
        // power
        mpoly2 = mpolynomial::MPolynomial::linear_pown(&[7.0, -3.0, 5.0], &[0, 1, 2], n_var, n);
        assert_eq!(mpoly1, mpoly2);
    }

    #[test]
    fn multiplication_as_power_no_const() {
        // compare (-3x1+5x2)^6 using MPolynomial::mult and MPolynomial::linear_pown
        let n_var = 2;
        let n = 6;
        // Create two copies of (3x1 +5x2)
        let mut mpoly1 = mpolynomial::MPolynomial::new(n_var);
        mpoly1.add(&vec![1, 0], -3.0);
        mpoly1.add(&vec![0, 1], 5.0);
        let mut mpoly2 = mpolynomial::MPolynomial::new(n_var);
        mpoly2.add(&vec![1, 0], -3.0);
        mpoly2.add(&vec![0, 1], 5.0);

        // multiply
        for _ in 1..n {
            mpoly1.mult(&mpoly2);
        }
        // power
        mpoly2 = mpolynomial::MPolynomial::linear_pown(&[-3.0, 5.0], &[1, 2], n_var, n);
        assert_eq!(mpoly1, mpoly2);
    }
    #[test]
    fn square() {
        // compare (-3x1+5x2)^6 using MPolynomial::mult and MPolynomial::linear_pown
        let n_var = 2;
        let n = 2;
        // Create two copies of (3x1 +5x2)
        let mut mpoly1 = mpolynomial::MPolynomial::new(n_var);
        mpoly1.add(&vec![1, 0], -3.0);
        mpoly1.add(&vec![0, 1], 5.0);
        mpoly1.square();

        // power
        let mpoly2 = mpolynomial::MPolynomial::linear_pown(&[-3.0, 5.0], &[1, 2], n_var, n);
        assert_eq!(mpoly1, mpoly2);
    }
    #[test]
    fn cube() {
        // compare (-3x1+5x2)^6 using MPolynomial::mult and MPolynomial::linear_pown
        let n_var = 2;
        let n = 3;
        // Create two copies of (3x1 +5x2)
        let mut mpoly1 = mpolynomial::MPolynomial::new(n_var);
        mpoly1.add(&vec![1, 0], -3.0);
        mpoly1.add(&vec![0, 1], 5.0);
        mpoly1.cube();

        // power
        let mpoly2 = mpolynomial::MPolynomial::linear_pown(&[-3.0, 5.0], &[1, 2], n_var, n);
        assert_eq!(mpoly1, mpoly2);
    }
    #[test]
    fn power() {
        // compare (-3x1+5x2)^6 using MPolynomial::mult and MPolynomial::linear_pown
        let n_var = 2;
        let n = 100;
        // Create two copies of (3x1 +5x2)
        let mut mpoly1 = mpolynomial::MPolynomial::new(n_var);
        mpoly1.add(&vec![1, 0], -3.0);
        mpoly1.add(&vec![0, 1], 5.0);
        let mut mpoly2 = mpoly1.clone();

        // multiply
        for _ in 1..n {
            mpoly1.mult(&mpoly2);
        }

        // power 1
        mpoly2.pown(n);
        // power 2
        //let mpoly2 = mpolynomial::MPolynomial::linear_pown(&[-3.0, 5.0], &[1, 2], n_var, n);
        assert_eq!(mpoly1, mpoly2);
    }
    #[test]
    fn power2() {
        // compare (-3x1+5x2)^6 using MPolynomial::mult and MPolynomial::linear_pown
        let n_var = 3;
        let n = 3;
        // Create two copies of (3x1 +5x2)
        let mut mpoly1 = mpolynomial::MPolynomial::new(n_var);
        mpoly1.add(&vec![1, 0,0], -3.0);
        mpoly1.add(&vec![0, 1,0], 5.0);
        mpoly1.add(&vec![0, 0,1], 7.0);
        let mut mpoly2 = mpoly1.clone();

        // multiply
        for _ in 1..n {
            mpoly1.mult(&mpoly2);
        }

        // power 1
        mpoly2.pown2(n);
        // power 2
        //let mpoly2 = mpolynomial::MPolynomial::linear_pown(&[-3.0, 5.0], &[1, 2], n_var, n);
        println!("{}\n{}", mpoly1, mpoly2);
        assert_eq!(mpoly1, mpoly2);
    }
}
