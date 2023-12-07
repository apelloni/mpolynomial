use crate::rational::MPolyRat;
use crate::MPolynomial;
use num::{BigInt, BigRational};
use pest::iterators::Pairs;
use pest::Parser;

#[derive(Parser)]
#[grammar = "grammar.pest"]
pub struct INIParser;

pub fn parse_mpolyrat(expr: &str, var_names: &[String]) -> MPolyRat {
    let parsed = INIParser::parse(Rule::calculation, expr)
        .expect("Unsuccessful Parse") // unwrap the parse result
        .next()
        .unwrap(); // get and unwrap the `file` rule; never fails
    let n_var = var_names.len();
    let mut parsed_polynomial = parsed_eval_rat(parsed.into_inner(), n_var, var_names);
    parsed_polynomial.drop_zeros();
    parsed_polynomial
}
pub fn parse_polynomial(expr: &str, var_names: &[String]) -> MPolynomial<BigRational> {
    let parsed = INIParser::parse(Rule::calculation, expr)
        .expect("Unsuccessful Parse") // unwrap the parse result
        .next()
        .unwrap(); // get and unwrap the `file` rule; never fails
    let n_var = var_names.len();
    let mut parsed_polynomial = parsed_eval_poly(parsed.into_inner(), n_var, var_names);
    parsed_polynomial.drop_zeros();
    parsed_polynomial
}

pub fn parsed_eval_poly(
    parser_rules: Pairs<Rule>,
    n_var: usize,
    var_names: &[String],
) -> MPolynomial<BigRational> {
    let mut block = MPolynomial::new(n_var);
    let mut out = MPolynomial::new(n_var);
    let mut res = MPolynomial::new(n_var);
    let mut pows = vec![0; n_var];
    let mut exponent: i32;
    let mut last_rule = Rule::add;
    let mut start_of_input = true;
    let mut found_var;
    let mut dpower = false;
    for parsed_rule in parser_rules.into_iter() {
        //println!("{:?} \"{}\"", parsed_rule.as_rule(), parsed_rule.as_str());
        //println!("   Last rule: {:?}", last_rule);
        match parsed_rule.as_rule() {
            Rule::num | Rule::var | Rule::expr => {
                // Process adjacent expressions as if there is an
                // implicit multiplicaiton in between of them
                if last_rule == Rule::num || last_rule == Rule::var || last_rule == Rule::expr {
                    last_rule = Rule::multiply;
                }

                match (last_rule, parsed_rule.as_rule()) {
                    // Process the information in the expression in terms of polynomial
                    (Rule::power, Rule::num) => {
                        // Process the power of an expression considering
                        // that one power its already stored in the block
                        // so we need to generate the remaining (power-1)
                        // contirbutions
                        exponent = parsed_rule
                            .clone()
                            .into_inner()
                            .as_str()
                            .parse::<i32>()
                            .expect("Exponent Parsing Error");
                        //assert!(exponent > 0, "Exponents must be positive!");
                        if dpower {
                            if exponent <= 0 {
                                res.pown((-exponent + 1) as usize);
                                last_rule = Rule::multiply;
                            }
                            if exponent > 0 {
                                res.pown((exponent - 1) as usize);
                                last_rule = Rule::divide;
                            }
                        } else {
                            if exponent <= 0 {
                                res.pown((-exponent + 1) as usize);
                                last_rule = Rule::divide;
                            }
                            if exponent > 0 {
                                res.pown((exponent - 1) as usize);
                                last_rule = Rule::multiply;
                            }
                        }
                    }
                    (_, Rule::num) => {
                        let coeff = BigRational::from(
                            parsed_rule
                                .clone()
                                .into_inner()
                                .as_str()
                                .parse::<BigInt>()
                                .expect("BigInt Parsing Fail"),
                        );
                        for p in pows.iter_mut() {
                            *p = 0;
                        }
                        res.clear();
                        res.add_coeff(&pows, coeff);
                        dpower = last_rule == Rule::divide;
                    }
                    (_, Rule::var) => {
                        found_var = false;
                        let coeff: BigRational = num::one();
                        for (vn, p) in pows.iter_mut().enumerate() {
                            if parsed_rule.clone().into_inner().as_str() == var_names[vn] {
                                *p = 1;
                                found_var = true;
                            } else {
                                *p = 0;
                            }
                        }
                        if !found_var {
                            panic!(
                                "Unknown variable {:?} found. [ vars={:?} ]",
                                parsed_rule.into_inner().as_str(),
                                var_names
                            );
                        }
                        res.clear();
                        res.add_coeff(&pows, coeff);
                        dpower = last_rule == Rule::divide;
                    }
                    (_, Rule::expr) => {
                        res = parsed_eval_poly(parsed_rule.clone().into_inner(), n_var, var_names);
                        dpower = last_rule == Rule::divide;
                    }
                    _ => panic!("Impossible!"),
                };
                // append operation to the block based on `last_rule`
                match last_rule {
                    Rule::add => {
                        block = res.clone();
                    }
                    Rule::subtract => {
                        block.clear();
                        block -= &res;
                    }
                    Rule::multiply => {
                        block.mult(&res);
                    }
                    Rule::divide => {
                        if res.is_constant() {
                            block.scale(&(num::one::<BigRational>() / res.coeffs[0].clone()));
                        } else {
                            block
                                .exact_division(&res)
                                .expect("Monomial does not factorize!");
                        }
                    }
                    _ => {
                        panic!("Unknown operation sequence! .. {:?}  expr ..", last_rule)
                    }
                };
                last_rule = Rule::expr;
            }
            // When we have addition or subtraction we evaluate the block and add it to
            // the output polynomial
            Rule::add | Rule::subtract => {
                match last_rule {
                    Rule::num | Rule::var | Rule::expr => {
                        if start_of_input {
                            start_of_input = false;
                            out = block.clone();
                        } else {
                            out += &block;
                        }
                        block.clear();
                        last_rule = parsed_rule.as_rule();
                    }
                    Rule::add => {
                        last_rule = if parsed_rule.as_rule() == Rule::add {
                            Rule::add
                        } else {
                            Rule::subtract
                        };
                    }
                    Rule::subtract => {
                        last_rule = if parsed_rule.as_rule() == Rule::subtract {
                            Rule::add
                        } else {
                            Rule::subtract
                        };
                    }
                    _ => {
                        panic!(
                            "Unknown operation sequence! .. {:?}  {:?} ..",
                            last_rule,
                            parsed_rule.as_rule()
                        )
                    }
                };
            }
            // The information about the binary operations *|/|^ that are performed within
            // the block are used only a posteriori as information in `last_rule`
            // The double multiplication is treated as a power operator
            Rule::multiply | Rule::divide | Rule::power => {
                match last_rule {
                    Rule::num | Rule::var | Rule::expr => {
                        last_rule = parsed_rule.as_rule();
                    }
                    Rule::multiply => last_rule = Rule::power,
                    _ => {
                        panic!(
                            "Unknown operation sequence! .. {:?}  {:?} ..",
                            last_rule,
                            parsed_rule.as_rule()
                        )
                    }
                };
            }
            x => {
                println!("{:?}", x);
                println!(" {:?}", parsed_rule.as_span());
                last_rule = x;
            }
        };
        //println!(
        //    " -> {} [{}]",
        //    out.to_str(var_names),
        //    block.to_str(var_names)
        //);
    }
    if block.coeffs.len() != 0 {
        out += &block;
    }
    //println!("{}", out);
    //println!("{}", out.to_str(var_names));
    out
}

pub fn parsed_eval_rat(parser_rules: Pairs<Rule>, n_var: usize, var_names: &[String]) -> MPolyRat {
    let mut block = MPolyRat::new(n_var);
    let mut out = MPolyRat::new(n_var);
    let mut res = MPolyRat::new(n_var);
    let mut pows = vec![0; n_var];
    let mut exponent: i32;
    let mut last_rule = Rule::add;
    let mut start_of_input = true;
    let mut found_var;
    let mut dpower = false;
    for parsed_rule in parser_rules.into_iter() {
        //println!("{:?} \"{}\"", parsed_rule.as_rule(), parsed_rule.as_str());
        match parsed_rule.as_rule() {
            Rule::num | Rule::var | Rule::expr => {
                // Process adjacent expressions as if there is an
                // implicit multiplicaiton in between of them
                if last_rule == Rule::num || last_rule == Rule::var || last_rule == Rule::expr {
                    last_rule = Rule::multiply;
                }

                match (last_rule, parsed_rule.as_rule()) {
                    // Process the information in the expression in terms of polynomial
                    (Rule::power, Rule::num | Rule::expr) => {
                        // Process the power of an expression considering
                        // that one power its already stored in the block
                        // so we need to generate the remaining (power-1)
                        // contirbutions
                        exponent = parsed_rule
                            .clone()
                            .into_inner()
                            .as_str()
                            .parse::<i32>()
                            .expect("Exponent Parsing Error");
                        //assert!(exponent > 0, "Exponents must be positive!");
                        if dpower {
                            if exponent <= 0 {
                                res.pown((-exponent + 1) as usize);
                                last_rule = Rule::multiply;
                            }
                            if exponent > 0 {
                                res.pown((exponent - 1) as usize);
                                last_rule = Rule::divide;
                            }
                        } else {
                            if exponent <= 0 {
                                res.pown((-exponent + 1) as usize);
                                last_rule = Rule::divide;
                            }
                            if exponent > 0 {
                                res.pown((exponent - 1) as usize);
                                last_rule = Rule::multiply;
                            }
                        }
                    }
                    (_, Rule::num) => {
                        let coeff = BigRational::from(
                            parsed_rule
                                .clone()
                                .into_inner()
                                .as_str()
                                .parse::<BigInt>()
                                .expect("BigInt Parsing Fail"),
                        );
                        for p in pows.iter_mut() {
                            *p = 0;
                        }
                        res.clear();
                        res.num.add_coeff(&pows, coeff);
                        dpower = last_rule == Rule::divide;
                    }
                    (_, Rule::var) => {
                        found_var = false;
                        let coeff: BigRational = num::one();
                        for (vn, p) in pows.iter_mut().enumerate() {
                            if parsed_rule.clone().into_inner().as_str() == var_names[vn] {
                                *p = 1;
                                found_var = true;
                            } else {
                                *p = 0;
                            }
                        }
                        if !found_var {
                            panic!(
                                "Unknown variable {:?} found. [ vars={:?} ]",
                                parsed_rule.into_inner().as_str(),
                                var_names
                            );
                        }
                        res.clear();
                        res.num.add_coeff(&pows, coeff);
                        dpower = last_rule == Rule::divide;
                    }
                    (_, Rule::expr) => {
                        res = parsed_eval_rat(parsed_rule.into_inner(), n_var, var_names);
                        dpower = last_rule == Rule::divide;
                    }
                    _ => panic!("Impossible!"),
                };
                // append operation to the block based on `last_rule`
                match last_rule {
                    Rule::add => {
                        block = res.clone();
                    }
                    Rule::subtract => {
                        block.clear();
                        block -= &res;
                    }
                    Rule::multiply => {
                        block *= &res;
                    }
                    Rule::divide => {
                        if res.is_constant() {
                            block.scale(&(num::one::<BigRational>() / res.num.coeffs[0].clone()));
                        } else {
                            block /= &res;
                        }
                    }
                    _ => {
                        panic!("Unknown operation sequence! .. {:?}  expr ..", last_rule)
                    }
                };
                last_rule = Rule::expr;
            }
            // When we have addition or subtraction we evaluate the block and add it to
            // the output polynomial
            Rule::add | Rule::subtract => {
                match last_rule {
                    Rule::num | Rule::var | Rule::expr => {
                        if start_of_input {
                            start_of_input = false;
                            out = block.clone();
                        } else {
                            out += &block;
                        }
                        block.clear();
                        last_rule = parsed_rule.as_rule();
                    }
                    Rule::add => {
                        last_rule = if parsed_rule.as_rule() == Rule::add {
                            Rule::add
                        } else {
                            Rule::subtract
                        };
                    }
                    Rule::subtract => {
                        last_rule = if parsed_rule.as_rule() == Rule::subtract {
                            Rule::add
                        } else {
                            Rule::subtract
                        };
                    }
                    _ => {
                        panic!(
                            "Unknown operation sequence! .. {:?}  {:?} ..",
                            last_rule,
                            parsed_rule.as_rule()
                        )
                    }
                };
            }
            // The information about the binary operations *|/|^ that are performed within
            // the block are used only a posteriori as information in `last_rule`
            // The double multiplication is treated as a power operator
            Rule::multiply | Rule::divide | Rule::power => {
                match last_rule {
                    Rule::num | Rule::var | Rule::expr => {
                        last_rule = parsed_rule.as_rule();
                    }
                    Rule::multiply => last_rule = Rule::power,
                    _ => {
                        panic!(
                            "Unknown operation sequence! .. {:?}  {:?} ..",
                            last_rule,
                            parsed_rule.as_rule()
                        )
                    }
                };
            }
            x => {
                println!("{:?}", x);
                println!(" {:?}", parsed_rule.as_span());
                last_rule = x;
            }
        };
        //println!(
        //    " -> {} [{}]",
        //    out.to_str(var_names),
        //    block.to_str(var_names)
        //);
    }
    if block.num.coeffs.len() != 0 {
        out += &block;
    }
    //println!("{}", out);
    //println!("{}", out.to_str(var_names));
    out
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn parsing_test_polynomial() {
        //let expr = INIParser::parse(Rule::calculation, "w+a/ll123 + c4^3")
        let var_names = vec![String::from("w"), String::from("a")];

        //let expr_str = "w^0 ";
        let expr_str = "10+a^3*(7+a-a^2)*(w+1)^0";
        println!("input: {}", expr_str);
        let expr_parsed = parse_polynomial(expr_str, var_names.as_slice());
        println!("output: {}", expr_parsed.to_str(var_names.as_slice()));
        assert_eq!(
            "+10+7*a^3+a^4-a^5",
            expr_parsed.to_str(var_names.as_slice())
        );
        // TEST
        let expr_str = "1/2*(1-w/4+w*(1/2-3/4))*4+w^2*(w^23*a^2*w^5)";
        println!("input: {}", expr_str);
        let expr_parsed = parse_polynomial(expr_str, var_names.as_slice());
        println!("output: {}", expr_parsed.to_str(var_names.as_slice()));
        assert_eq!("+2-w^1+w^30*a^2", expr_parsed.to_str(var_names.as_slice()));

        let expr_str = "1/(1-45^3-4*w*w+(2*w)^2)";
        println!("input: {}", expr_str);
        let expr_parsed = parse_polynomial(expr_str, var_names.as_slice());
        println!("output: {}", expr_parsed.to_str(var_names.as_slice()));
        assert_eq!("-1/91124", expr_parsed.to_str(var_names.as_slice()));

        let expr_str = "0+(w)*(-1)+(-w)*(-1)+(-w)+(w)+(-2)*(-2)+(-w)+(2 + w)";
        println!("input: {}", expr_str);
        let expr_parsed = parse_polynomial(expr_str, var_names.as_slice());
        println!("output: {}", expr_parsed.to_str(var_names.as_slice()));
        assert_eq!("+6", expr_parsed.to_str(var_names.as_slice()));

        let expr_str = "w^0";
        println!("input: {}", expr_str);
        let expr_parsed = parse_polynomial(expr_str, var_names.as_slice());
        println!("output: {}", expr_parsed.to_str(var_names.as_slice()));
        assert_eq!("+1", expr_parsed.to_str(var_names.as_slice()));
    }

    #[test]
    fn parsing_test_mpolyrat() {
        //let expr = INIParser::parse(Rule::calculation, "w+a/ll123 + c4^3")
        let var_names = vec![String::from("w"), String::from("a")];

        let expr_str = "(1-w)^2/(1-w)";
        let mut expr_parsed = parse_mpolyrat(expr_str, var_names.as_slice());
        assert_eq!(
            "(+1-2*w^1+w^2)/(+1-w^1)",
            expr_parsed.to_str(var_names.as_slice())
        );
        expr_parsed.reduce();
        assert_eq!("(+1-w^1)/(+1)", expr_parsed.to_str(var_names.as_slice()));

        //let expr_str = "w^0 ";
        let expr_str = "10+a^3*(7+a-a^2)*(w+1)^0";
        println!("input: {}", expr_str);
        let mut expr_parsed = parse_mpolyrat(expr_str, var_names.as_slice());
        expr_parsed.reduce();
        println!("output: {}", expr_parsed.to_str(var_names.as_slice()));
        assert_eq!(
            "(+10+7*a^3+a^4-a^5)/(+1)",
            expr_parsed.to_str(var_names.as_slice())
        );
        // TEST
        let expr_str = "1/2*(1-w/4+w*(1/2-3/4))*4+w^2*(w^23*a^2*w^5)";
        println!("input: {}", expr_str);
        let expr_parsed = parse_mpolyrat(expr_str, var_names.as_slice());
        assert_eq!(
            "(+2-w^1+w^30*a^2)/(+1)",
            expr_parsed.to_str(var_names.as_slice())
        );

        let expr_str = "1/(1-45^3-4*w*w+(2*w)^2)";
        println!("input: {}", expr_str);
        let expr_parsed = parse_mpolyrat(expr_str, var_names.as_slice());
        println!("output: {}", expr_parsed.to_str(var_names.as_slice()));
        assert_eq!("(-1/91124)/(+1)", expr_parsed.to_str(var_names.as_slice()));

        let expr_str = "0+(w)*(-1)+(-w)*(-1)+(-w)+(w)+(-2)*(-2)+(-w)+(2 + w)";
        println!("input: {}", expr_str);
        let expr_parsed = parse_mpolyrat(expr_str, var_names.as_slice());
        println!("output: {}", expr_parsed.to_str(var_names.as_slice()));
        assert_eq!("(+6)/(+1)", expr_parsed.to_str(var_names.as_slice()));

        let expr_str = "w^0";
        println!("input: {}", expr_str);
        let mut expr_parsed = parse_mpolyrat(expr_str, var_names.as_slice());
        expr_parsed.reduce();
        println!("output: {}", expr_parsed.to_str(var_names.as_slice()));
    }

    #[test]
    fn minus_sign() {
        //let expr = INIParser::parse(Rule::calculation, "w+a/ll123 + c4^3")
        let var_names = vec![String::from("w"), String::from("a")];

        // Check spaces
        let expr_str = "-(1-w)";
        let expr_parsed = parse_polynomial(expr_str, var_names.as_slice());
        assert_eq!("-1+w^1", expr_parsed.to_str(var_names.as_slice()));
    }

    #[test]
    fn whitespaces() {
        //let expr = INIParser::parse(Rule::calculation, "w+a/ll123 + c4^3")
        let var_names = vec![String::from("w"), String::from("a")];

        // Check spaces
        let expr_str = "( -     2  *           w  ^  3  + 4 )";
        let expr_parsed = parse_polynomial(expr_str, var_names.as_slice());
        assert_eq!("+4-2*w^3", expr_parsed.to_str(var_names.as_slice()));
    }

    #[test]
    fn implicit_multiplication() {
        //let expr = INIParser::parse(Rule::calculation, "w+a/ll123 + c4^3")
        let var_names = vec![String::from("w"), String::from("a")];

        // Check spaces
        let expr_str = "10a^3*2(7+a-a^2)*(w+1)^0";
        let expr_parsed = parse_polynomial(expr_str, var_names.as_slice());
        assert_eq!(
            "+140*a^3+20*a^4-20*a^5",
            expr_parsed.to_str(var_names.as_slice())
        );
        // TEST
        let expr_str = "1/2(1-w/4+w(1/2-3/4))4+w^2*(w^23a^2w^5)";
        let expr_parsed = parse_polynomial(expr_str, var_names.as_slice());
        assert_eq!("+2-w^1+w^30*a^2", expr_parsed.to_str(var_names.as_slice()));

        let expr_str = "( -     2             w  ^  3  + 4 )";
        let expr_parsed = parse_polynomial(expr_str, var_names.as_slice());
        assert_eq!("+4-2*w^3", expr_parsed.to_str(var_names.as_slice()));
    }

    #[test]
    fn power_as_double_star() {
        //let expr = INIParser::parse(Rule::calculation, "w+a/ll123 + c4^3")
        let var_names = vec![String::from("w"), String::from("a")];

        let expr_str = "10+a**3*(7+a-a**2)*(w+1)**0";
        let expr_parsed = parse_polynomial(expr_str, var_names.as_slice());
        //println!("input: {}", expr_str);
        //println!("output: {}", expr_parsed.to_str(var_names.as_slice()));
        assert_eq!(
            "+10+7*a^3+a^4-a^5",
            expr_parsed.to_str(var_names.as_slice())
        );

        let expr_str = "1/2*(1-w/4+w*(1/2-3/4))*4+w**2*(w**23*a**2*w**5)";
        let expr_parsed = parse_polynomial(expr_str, var_names.as_slice());
        assert_eq!("+2-w^1+w^30*a^2", expr_parsed.to_str(var_names.as_slice()));

        let expr_str = "1/(1-45**3-4*w*w+(2*w)**2)";
        let expr_parsed = parse_polynomial(expr_str, var_names.as_slice());
        assert_eq!("-1/91124", expr_parsed.to_str(var_names.as_slice()));

        let expr_str = "w**0";
        let expr_parsed = parse_polynomial(expr_str, var_names.as_slice());
        assert_eq!("+1", expr_parsed.to_str(var_names.as_slice()));
    }

    #[test]
    fn parse_power_polynomial_num() {
        let var_names = vec![String::from("x1"), String::from("x2")];
        let expr_str = "x1/2^3";
        let res_str = "1/8*x1";
        let expr = parse_polynomial(expr_str, &var_names);
        let res = parse_polynomial(res_str, &var_names);
        println!(">> {}", expr);
        assert_eq!(res, expr);
    }

    #[test]
    fn parse_dpower_polynomial_var() {
        let var_names = vec![String::from("x1"), String::from("x2")];
        let expr_str = "x1^4/x1^3";
        let res_str = "x1";
        let expr = parse_polynomial(expr_str, &var_names);
        let res = parse_polynomial(res_str, &var_names);
        println!(">> {}", expr);
        assert_eq!(res, expr);
    }

    #[test]
    fn parse_dpower_polynomial_expr() {
        let var_names = vec![String::from("x1"), String::from("x2")];
        let expr_str = "x1/(2+1)^3";
        let res_str = "1/27*x1";
        let expr = parse_polynomial(expr_str, &var_names);
        let res = parse_polynomial(res_str, &var_names);
        println!(">> {}", expr);
        assert_eq!(res, expr);
    }

    #[test]
    fn parse_dpower_polynomial_sum() {
        let var_names = vec![String::from("x1"), String::from("x2")];
        let expr_str = "x1/(2+1)^3+x1*9^2/3/3+x2^4";
        let res_str = "1/27*x1+9*x1+x2^4";
        let expr = parse_polynomial(expr_str, &var_names);
        let res = parse_polynomial(res_str, &var_names);
        println!(">> {}", expr);
        println!(">> {}", res);
        assert_eq!(res, expr);
    }

    #[test]
    fn parse_dpower_mpolyrat_num() {
        let var_names = vec![String::from("x1"), String::from("x2")];
        let expr_str = "(1-x1)/(1+3)^3";
        let res_str = "1/4/4/4*(1-x1)";
        let expr = parse_mpolyrat(expr_str, &var_names);
        let res = parse_mpolyrat(res_str, &var_names);
        println!(">> {}", expr);
        println!(">> {}", res);
        assert_eq!(res, expr);
    }

    #[test]
    fn parse_dpower_mpolyrat_var() {
        let var_names = vec![String::from("x1"), String::from("x2")];
        let expr_str = "(1-x1)/(x2)^2";
        let res_str = "(1-x1)/x2/x2";
        let expr = parse_mpolyrat(expr_str, &var_names);
        let res = parse_mpolyrat(res_str, &var_names);
        println!(">> {}", expr);
        println!(">> {}", res);
        assert_eq!(res, expr);
    }

    #[test]
    fn parse_dpower_mpolyrat_expr() {
        let var_names = vec![String::from("x1"), String::from("x2")];
        let expr_str = "(1-x1)^2/(x2+1)^2*(x2-1)^2";
        let res_str = "(1-x1)*(1-x1)/(x2+1)/(x2+1)*(x2-1)*(x2-1)";
        let expr = parse_mpolyrat(expr_str, &var_names);
        let res = parse_mpolyrat(res_str, &var_names);
        println!(">> {}", expr);
        println!(">> {}", res);
        assert_eq!(res, expr);
    }

    #[test]
    fn geometric_division() {
        let var_names = vec![String::from("x1")];

        let mpoly_string = "(1-x1)^2/(1-x1)";
        let mpolyrat = parse_mpolyrat(mpoly_string, &var_names);
        println!(">> {}", mpolyrat);

        // Geometric Series
        let  mpoly_string = "(1-x1^11)^2/(1-x1)^2";
        let res_str = "(1+x1+x1^2+x1^3+x1^4+x1^5+x1^6+x1^7+x1^8+x1^9+x1^10)^2";
        let mut mpolyrat = parse_mpolyrat(mpoly_string, &var_names);
        mpolyrat.reduce();
        let res = parse_mpolyrat(res_str, &var_names);
        assert_eq!(res.to_str(&var_names), mpolyrat.to_str(&var_names));
    }

    #[test]
    fn exact_division() {
        let var_names = vec![
            String::from("x1"),
            String::from("x2"),
            String::from("x3"),
            String::from("x4"),
        ];

        // Geometric Series
        let mpoly_str = "(1-x1^11)^2";
        let factor_str = "(1-x1)^2";
        let res_str = "(1+x1+x1^2+x1^3+x1^4+x1^5+x1^6+x1^7+x1^8+x1^9+x1^10)^2";

        let mut mpoly = parse_polynomial(mpoly_str, &var_names);
        let factor = parse_polynomial(factor_str, &var_names);
        let res = parse_polynomial(res_str, &var_names);

        mpoly.exact_division(&factor).unwrap();
        println!(">> {}", mpoly);
        assert_eq!(res.to_str(&var_names), mpoly.to_str(&var_names));

        // Multivariate
        let mpoly_str = "(11 x1^32 + 12 x2 + 99 x3 + 21 x1*x4)^10";
        let factor_str = "(11 x1**32 + 12 x2 + 99 x3 + 21 x1*x4)^8";

        let factor = parse_polynomial(factor_str, &var_names);
        let mut mpoly = parse_polynomial(mpoly_str, &var_names);

        mpoly.exact_division(&factor).unwrap();
        println!(">> {}", mpoly);
        let mut ss = String::new();
        ss += "+9801*x3^2+2376*x2^1*x3^1+144*x2^2+4158*x1^1*x3^1*x4^1";
        ss += "+504*x1^1*x2^1*x4^1+441*x1^2*x4^2+2178*x1^32*x3^1";
        ss += "+264*x1^32*x2^1+462*x1^33*x4^1+121*x1^64";
        assert_eq!(ss, mpoly.to_str(&var_names));
    }

    #[test]
    fn univariate_rational() {
        let variables = &[
            "x1".to_string(),
            "x2".to_string(),
            "x3".to_string(),
            "x4".to_string(),
        ];

        let mut mpolyrat = parse_mpolyrat(
            "27 / (29 * (37 / (28 * x1) + (28 * x1) / 2479)) + 37 / (28 * x1) + (28 * x1) / 2479",
            variables,
        );
        println!("{}", mpolyrat.to_str(variables));
        mpolyrat.reduce();
        println!("{}", mpolyrat.to_str(variables));
        assert_eq!(
            "(+37/28+52276/71891*x1^2+21952/227381317*x1^4)/(+x1^1+784/91723*x1^3)",
            mpolyrat.to_str(variables)
        );
    }

    #[test]
    fn multivariate_rational() {
        let variables = &[
            "x1".to_string(),
            "x2".to_string(),
            "x3".to_string(),
            "x4".to_string(),
            "a".to_string(),
        ];

        let mut input_string = String::new();
        input_string += "11905/(1904*a)+(27368*a)/9849+x1+x2+x3^(-1)+2/((2*a)/65+x1^(-1)";
        input_string += "+x2^(-1)+x3^(-1)+x4^(-1))+((77*a)/16+x1^(-1)+x2+x3+x4^(-1))^";
        input_string += "(-1)+x4^(-1)+2*((10*a)/43+x1^(-1)+x2+x3^(-1)+x4)+((35*a)/34+";
        input_string += "x1+x2^(-1)+x3+x4)^(-1)+2*((2*a)/7+x1^(-1)+x2+x3+x4)+((79*";
        input_string += "a)/30+x1^(-1)+x2+x3+x4)^(-1)+2*(25/(46*a)+x1+x2+x3+x4)+";
        input_string += "(99/(28*a)+(150603*a)/78352+2/x1+x1+3*x2+2/x3+x3+x4^(-1)+2*";
        input_string += "x4+(27/(89*a)+x1^(-1)+x2^(-1)+x3^(-1)+x4)^(-1))^(-1)+(23231/(11";
        input_string += "780*a)+(91*a)/24+2/(47/(24*a)+x1^(-1)+x2^(-1)+x3^(-1)+x4^(-1))";
        input_string += "+(31/(4*a)+x1+x2^(-1)+x3^(-1)+x4^(-1))^(-1)+2/((100*a)/47+x1^(";
        input_string += "-1)+x2+x3+x4^(-1))+2*((9*a)/7+x1+x2+x3+x4^(-1))+3*(9/(4";
        input_string += "7*a)+x1^(-1)+x2^(-1)+x3^(-1)+x4)+2*((31*a)/94+x1^(-1)+x2+x3";
        input_string += "^(-1)+x4)+((27*a)/2+x1^(-1)+x2^(-1)+x3+x4)^(-1)+3*((34*a)/71";
        input_string += "+x1+x2^(-1)+x3+x4))^(-1)+(4109/(680*a)+(10*a)/13+2*x1+x2^(-1";
        input_string += ")+x2+2/x3+2/(31/(25*a)+x1^(-1)+x2+x3^(-1)+x4^(-1))+((16*a)/";
        input_string += "45+x1^(-1)+x2+x3+x4^(-1))^(-1)+x4^(-1)+x4+2/(a/5+x1+x2^(-";
        input_string += "1)+x3^(-1)+x4)+((63*a)/62+x1+x2+x3+x4)^(-1)+(2868539/(604";
        input_string += "219*a)+x1^(-1)+2*x1+x2^(-1)+2*x2+3/x3+x4^(-1)+2*x4+(49/(39*";
        input_string += "a)+x1^(-1)+x2^(-1)+x3+x4)^(-1))^(-1)+(335/(24*a)+(17*a)/21+";
        input_string += "x1^(-1)+x1+x2^(-1)+x2+x3^(-1)+x3+x4^(-1)+x4+((19*a)/3+x1^";
        input_string += "(-1)+x2+x3^(-1)+x4)^(-1)+(71/(64*a)+x1^(-1)+x2+x3+x4)^(-1))";
        input_string += "^(-1))^(-1)+(195047/(38335*a)+(3367*a)/2788+x1^(-1)+x2+x3+2*((3";
        input_string += "5*a)/69+x1^(-1)+x2^(-1)+x3^(-1)+x4^(-1))+2*(64/(93*a)+x1+x2^(";
        input_string += "-1)+x3^(-1)+x4^(-1))+(31/(22*a)+x1^(-1)+x2^(-1)+x3+x4^(-1))^(";
        input_string += "-1)+x4^(-1)+((9*a)/68+x1^(-1)+x2^(-1)+x3^(-1)+x4)^(-1)+(4/(";
        input_string += "77*a)+x1+x2+x3^(-1)+x4)^(-1)+2*(100/(83*a)+x1^(-1)+x2+x3+";
        input_string += "x4)+2/(54/(85*a)+x1+x2+x3+x4)+((23*a)/4+(27/(8*a)+x1^(-1)";
        input_string += "+x2^(-1)+x3^(-1)+x4^(-1))^(-1)+((43*a)/75+x1+x2+x3+x4^(-1))";
        input_string += "^(-1)+(12*a+x1+x2+x3^(-1)+x4)^(-1)+(7*a+x1+x2+x3+x4)^(-";
        input_string += "1))^(-1))^(-1)27/(29*(37/(28*x1)+(28*x1)/2479))+37/(28*";
        input_string += "x1)+(28*x1)/2479";
        let mut mpolyrat = parse_mpolyrat(&input_string, variables);
        println!("{}", mpolyrat.to_str(variables));
        mpolyrat.reduce();
        println!("{}", mpolyrat.to_str(variables));
        assert_eq!(
            "(+37/28+52276/71891*x1^2+21952/227381317*x1^4)/(+x1^1+784/91723*x1^3)",
            mpolyrat.to_str(variables)
        );
    }

    #[test]
    fn overflow() {
        let var_names = vec![String::from("x1"), String::from("x2")];

        let mpoly_str = "(x1 + x2 )^200";
        println!("Check overflow for : {}", mpoly_str);
        let mut mpoly_a = MPolynomial::linear_pown(
            &[num::one::<BigRational>(), num::one::<BigRational>()],
            &[1],
            2,
            200,
        );
        mpoly_a.replace(
            1,
            &[num::one::<BigRational>(), num::one::<BigRational>()],
            &[1, 2],
        );

        let mpoly_b = parse_polynomial(mpoly_str, &var_names);
        assert_eq!(mpoly_b.to_str(&var_names), mpoly_a.to_str(&var_names));
    }
}
