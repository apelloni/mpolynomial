use crate::MPolynomial;
use num::{BigInt, BigRational};
use pest::iterators::Pairs;
use pest::Parser;

#[derive(Parser)]
#[grammar = "grammar.pest"]
pub struct INIParser;

pub fn parse_expression(expr: &str, var_names: &[String]) -> MPolynomial<BigRational> {
    let parsed = INIParser::parse(Rule::calculation, expr)
        .expect("Unsuccessful Parse") // unwrap the parse result
        .next()
        .unwrap(); // get and unwrap the `file` rule; never fails
    let n_var = var_names.len();
    let mut parsed_polynomial = parsed_eval_fun(parsed.into_inner(), n_var, var_names);
    parsed_polynomial.drop_zeros();
    parsed_polynomial
}

pub fn parsed_eval_fun(
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
                        if exponent <= 0 {
                            res.pown((-exponent + 1) as usize);
                            last_rule = Rule::divide;
                        }
                        if exponent > 0 {
                            res.pown((exponent - 1) as usize);
                            last_rule = Rule::multiply;
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
                        res.add(&pows, coeff);
                    }
                    (_, Rule::var) => {
                        let coeff: BigRational = num::one();
                        for (vn, p) in pows.iter_mut().enumerate() {
                            if parsed_rule.clone().into_inner().as_str() == var_names[vn] {
                                *p = 1;
                            } else {
                                *p = 0;
                            }
                        }
                        res.clear();
                        res.add(&pows, coeff);
                    }
                    (_, Rule::expr) => {
                        res = parsed_eval_fun(parsed_rule.into_inner(), n_var, var_names);
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

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn parsing_test() {
        //let expr = INIParser::parse(Rule::calculation, "w+a/ll123 + c4^3")
        let var_names = vec![String::from("w"), String::from("a")];

        //let expr_str = "w^0 ";
        let expr_str = "10+a^3*(7+a-a^2)*(w+1)^0";
        println!("input: {}", expr_str);
        let expr_parsed = parse_expression(expr_str, var_names.as_slice());
        println!("output: {}", expr_parsed.to_str(var_names.as_slice()));
        assert_eq!(
            "+10+7*a^3+a^4-a^5",
            expr_parsed.to_str(var_names.as_slice())
        );
        // TEST
        let expr_str = "1/2*(1-w/4+w*(1/2-3/4))*4+w^2*(w^23*a^2*w^5)";
        println!("input: {}", expr_str);
        let expr_parsed = parse_expression(expr_str, var_names.as_slice());
        println!("output: {}", expr_parsed.to_str(var_names.as_slice()));
        assert_eq!("+2-w^1+w^30*a^2", expr_parsed.to_str(var_names.as_slice()));

        let expr_str = "1/(1-45^3-4*w*w+(2*w)^2)";
        println!("input: {}", expr_str);
        let expr_parsed = parse_expression(expr_str, var_names.as_slice());
        println!("output: {}", expr_parsed.to_str(var_names.as_slice()));
        assert_eq!("-1/91124", expr_parsed.to_str(var_names.as_slice()));

        let expr_str = "0+(w)*(-1)+(-w)*(-1)+(-w)+(w)+(-2)*(-2)+(-w)+(2 + w)";
        println!("input: {}", expr_str);
        let expr_parsed = parse_expression(expr_str, var_names.as_slice());
        println!("output: {}", expr_parsed.to_str(var_names.as_slice()));
        assert_eq!("+6", expr_parsed.to_str(var_names.as_slice()));

        let expr_str = "w^0";
        println!("input: {}", expr_str);
        let expr_parsed = parse_expression(expr_str, var_names.as_slice());
        println!("output: {}", expr_parsed.to_str(var_names.as_slice()));
        assert_eq!("+1", expr_parsed.to_str(var_names.as_slice()));
    }

    #[test]
    fn minus_sign() {
        //let expr = INIParser::parse(Rule::calculation, "w+a/ll123 + c4^3")
        let var_names = vec![String::from("w"), String::from("a")];

        // Check spaces
        let expr_str = "-(1-w)";
        let expr_parsed = parse_expression(expr_str, var_names.as_slice());
        assert_eq!("-1+w^1", expr_parsed.to_str(var_names.as_slice()));
    }

    #[test]
    fn whitespaces() {
        //let expr = INIParser::parse(Rule::calculation, "w+a/ll123 + c4^3")
        let var_names = vec![String::from("w"), String::from("a")];

        // Check spaces
        let expr_str = "( -     2  *           w  ^  3  + 4 )";
        let expr_parsed = parse_expression(expr_str, var_names.as_slice());
        assert_eq!("+4-2*w^3", expr_parsed.to_str(var_names.as_slice()));
    }

    #[test]
    fn implicit_multiplication() {
        //let expr = INIParser::parse(Rule::calculation, "w+a/ll123 + c4^3")
        let var_names = vec![String::from("w"), String::from("a")];

        // Check spaces
        let expr_str = "10a^3*2(7+a-a^2)*(w+1)^0";
        let expr_parsed = parse_expression(expr_str, var_names.as_slice());
        assert_eq!(
            "+140*a^3+20*a^4-20*a^5",
            expr_parsed.to_str(var_names.as_slice())
        );
        // TEST
        let expr_str = "1/2(1-w/4+w(1/2-3/4))4+w^2*(w^23a^2w^5)";
        let expr_parsed = parse_expression(expr_str, var_names.as_slice());
        assert_eq!("+2-w^1+w^30*a^2", expr_parsed.to_str(var_names.as_slice()));

        let expr_str = "( -     2             w  ^  3  + 4 )";
        let expr_parsed = parse_expression(expr_str, var_names.as_slice());
        assert_eq!("+4-2*w^3", expr_parsed.to_str(var_names.as_slice()));
    }

    #[test]
    fn power_as_double_star() {
        //let expr = INIParser::parse(Rule::calculation, "w+a/ll123 + c4^3")
        let var_names = vec![String::from("w"), String::from("a")];

        let expr_str = "10+a**3*(7+a-a**2)*(w+1)**0";
        let expr_parsed = parse_expression(expr_str, var_names.as_slice());
        //println!("input: {}", expr_str);
        //println!("output: {}", expr_parsed.to_str(var_names.as_slice()));
        assert_eq!(
            "+10+7*a^3+a^4-a^5",
            expr_parsed.to_str(var_names.as_slice())
        );

        let expr_str = "1/2*(1-w/4+w*(1/2-3/4))*4+w**2*(w**23*a**2*w**5)";
        let expr_parsed = parse_expression(expr_str, var_names.as_slice());
        assert_eq!("+2-w^1+w^30*a^2", expr_parsed.to_str(var_names.as_slice()));

        let expr_str = "1/(1-45**3-4*w*w+(2*w)**2)";
        let expr_parsed = parse_expression(expr_str, var_names.as_slice());
        assert_eq!("-1/91124", expr_parsed.to_str(var_names.as_slice()));

        let expr_str = "w**0";
        let expr_parsed = parse_expression(expr_str, var_names.as_slice());
        assert_eq!("+1", expr_parsed.to_str(var_names.as_slice()));
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

        let mut mpoly = parse_expression(mpoly_str, &var_names);
        let factor = parse_expression(factor_str, &var_names);
        let res = parse_expression(res_str, &var_names);

        mpoly.exact_division(&factor).unwrap();
        assert_eq!(res.to_str(&var_names), mpoly.to_str(&var_names));

        // Multivariate
        let mpoly_str = "(11 x1^32 + 12 x2 + 99 x3 + 21 x1*x4)^10";
        let factor_str = "(11 x1**32 + 12 x2 + 99 x3 + 21 x1*x4)^8";

        let factor = parse_expression(factor_str, &var_names);
        let mut mpoly = parse_expression(mpoly_str, &var_names);

        mpoly.exact_division(&factor).unwrap();
        //println!("{}", mpoly);
        let mut ss = String::new();
        ss += "+9801*x3^2+2376*x2^1*x3^1+144*x2^2+4158*x1^1*x3^1*x4^1";
        ss += "+504*x1^1*x2^1*x4^1+441*x1^2*x4^2+2178*x1^32*x3^1";
        ss += "+264*x1^32*x2^1+462*x1^33*x4^1+121*x1^64";
        assert_eq!(ss, format!("{}", mpoly));
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

        let mpoly_b = parse_expression(mpoly_str, &var_names);
        assert_eq!(mpoly_b.to_str(&var_names), mpoly_a.to_str(&var_names));
    }
}
