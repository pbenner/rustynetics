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

use futures::executor::block_on;
use futures::Stream;
use futures::stream::StreamExt;

/* -------------------------------------------------------------------------- */

pub fn stream_to_blocking_iter<S>(mut stream: S) -> impl Iterator<Item = S::Item>
where
    S: Stream + Unpin, // Unpin is needed because we are passing the stream by mutable reference
{
    std::iter::from_fn(move || block_on(stream.next()))
}
