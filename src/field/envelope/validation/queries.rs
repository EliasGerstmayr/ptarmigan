//! Strict reader for the unquoted numeric CSV emitted by Report::queries.
use super::{gaussian::Gaussian, grid, runner::Query, Result};
use crate::field::envelope::UniformAxis;

pub fn parse(text: &str, g: &Gaussian, cells: &[usize], extent: f64) -> Result<Vec<Query>> {
    let mut lines = text.lines();
    let columns: Vec<_> = lines.next().ok_or("empty query CSV")?.split(',').collect();
    let mut indices = [0; 5];
    for (i, name) in ["x_m", "y_m", "z_m", "xi_m", "ct_m"].iter().enumerate() {
        let matches: Vec<_> = columns
            .iter()
            .enumerate()
            .filter(|(_, c)| *c == name)
            .collect();
        if matches.len() != 1 {
            return Err(format!("query CSV requires exactly one {} column", name).into());
        }
        indices[i] = matches[0].0;
    }
    let mut all_axes = Vec::new();
    for &n in cells {
        let nodes = grid::axes(g, n, extent)?;
        all_axes.push(
            nodes
                .iter()
                .map(|a| UniformAxis::try_from_coordinates(a))
                .collect::<std::result::Result<Vec<_>, _>>()?,
        );
    }
    let mut result = Vec::new();
    for (line, row) in lines.enumerate() {
        let number = line + 2;
        let fields: Vec<_> = row.split(',').collect();
        if fields.len() != columns.len() {
            return Err(format!("query CSV line {}: wrong column count", number).into());
        }
        let mut v = [0.0; 5];
        for i in 0..5 {
            v[i] = fields[indices[i]]
                .parse::<f64>()
                .map_err(|_| format!("query CSV line {}: invalid number", number))?;
        }
        if v.iter().any(|x| !x.is_finite()) {
            return Err(format!("query CSV line {}: nonfinite coordinate", number).into());
        }
        // Report uses 17 digits after the decimal in scientific notation. This
        // round-trips f64, including the effective xi previously computed as ct-z.
        // Require exact numerical equality (signed zeros compare equal), never
        // move points or accept a tolerance that could cross canonical bounds.
        if v[4] - v[2] != v[3] {
            return Err(format!(
                "query CSV line {}: xi must equal ct-z exactly after f64 parsing",
                number
            )
            .into());
        }
        for (resolution, axes) in cells.iter().zip(&all_axes) {
            for a in 0..4 {
                if v[a] < axes[a].coordinate(0).unwrap()
                    || v[a] > axes[a].coordinate(axes[a].len() - 1).unwrap()
                {
                    return Err(format!(
                        "query CSV line {}: axis {} outside canonical domain at {} cells",
                        number, a, resolution
                    )
                    .into());
                }
            }
        }
        result.push(Query {
            stored: [v[0], v[1], v[2], v[3]],
            lab: [v[4], v[0], v[1], v[2]],
        });
    }
    if result.is_empty() {
        return Err("query CSV contains no points".into());
    }
    Ok(result)
}

#[test]
fn saved_queries_roundtrip_and_reject_invalid_data() {
    use super::{gaussian::Parameters, runner::queries};
    use crate::field::Polarization;
    let g = Gaussian::new(Parameters::default(), Polarization::Linear).unwrap();
    let points = queries(&g, &[2, 3, 4], 2.0, 8, 17).unwrap();
    let header = "x_m,y_m,z_m,xi_m,ct_m\n";
    let mut csv = header.to_owned();
    for p in &points {
        csv += &format!(
            "{:.17e},{:.17e},{:.17e},{:.17e},{:.17e}\n",
            p.stored[0], p.stored[1], p.stored[2], p.stored[3], p.lab[0]
        );
    }
    assert_eq!(parse(&csv, &g, &[2, 3, 4], 2.0).unwrap(), points);
    for bad in [
        "x_m,y_m,z_m,ct_m\n0,0,0,0\n".to_owned(),
        header.to_owned(),
        format!("{}NaN,0,0,0,0\n", header),
        format!("{}0,0,0,1e-6,0\n", header),
        format!("{}1,0,0,0,0\n", header),
        format!("{}0,0,0,0\n", header),
        format!("{}0,0,0,inf,inf\n", header),
    ] {
        assert!(
            parse(&bad, &g, &[2, 3, 4], 2.0).is_err(),
            "accepted {}",
            bad
        );
    }
    // Choose a rounded domain whose canonical upper end depends on resolution.
    // A point admitted by the first grid must still fail a later grid's bounds.
    let mut found = false;
    for extent in [0.1, 0.3, 0.7, 1.1, 1.3, 1.7, 2.1] {
        let bounds: Vec<_> = (2..10)
            .map(|n| {
                let axes = grid::axes(&g, n, extent).unwrap();
                let a = UniformAxis::try_from_coordinates(&axes[0]).unwrap();
                (n, a.coordinate(a.len() - 1).unwrap())
            })
            .collect();
        if let Some(&(n, hi)) = bounds
            .iter()
            .find(|(_, hi)| bounds.iter().any(|(_, lo)| lo < hi))
        {
            let &(m, _) = bounds.iter().find(|(_, lo)| *lo < hi).unwrap();
            let csv = format!("{}{:.17e},0,0,0,0\n", header, hi);
            assert!(parse(&csv, &g, &[n], extent).is_ok());
            assert!(parse(&csv, &g, &[n, m], extent).is_err());
            found = true;
            break;
        }
    }
    assert!(found, "fixture must exercise different canonical endpoints");
}
