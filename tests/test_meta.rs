/* Copyright (C) 2024 Philipp Benner
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/* -------------------------------------------------------------------------- */

#[cfg(test)]
mod tests {

    use rustynetics::meta::{Meta};

    #[test]
    fn test_meta() {

        let names = vec!["name".to_string(), "age".to_string()];
        let data: Vec<Box<dyn MetaData>> = vec![
            Box::new(vec!["Alice".to_string(), "Bob".to_string()]),
            Box::new(vec![25, 30]),
        ];
        let meta = Meta::new(names, data).unwrap();
        println!("{}", meta);

    }
}
