/*

 Ripser: a lean C++ code for computation of Vietoris-Rips persistence barcodes

 MIT License

 Copyright (c) 2015–2019 Ulrich Bauer

 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:

 The above copyright notice and this permission notice shall be included in all
 copies or substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 SOFTWARE.

 You are under no obligation whatsoever to provide any bug fixes, patches, or
 upgrades to the features, functionality or performance of the source code
 ("Enhancements") to anyone; however, if you choose to make your Enhancements
 available either publicly, or directly to the author of this software, without
 imposing a separate written license agreement for such Enhancements, then you
 hereby grant the following license: a non-exclusive, royalty-free perpetual
 license to install, use, modify, prepare derivative works, incorporate into
 other computer software, distribute, and sublicense such enhancements or
 derivative works thereof, in binary and source code form.

*/

#include <benchmark/benchmark.h>

#define MAX_NUMBER_THREADS (48)

#include <chrono>
using std::chrono::duration;
using std::chrono::duration_cast;
using std::chrono::high_resolution_clock;
using std::chrono::milliseconds;

#include "../gph/src/ripser.h"

struct euclidean_distance_matrix {
    std::vector<std::vector<value_t>> points;
    std::vector<value_t> diagonal;

    euclidean_distance_matrix(std::vector<std::vector<value_t>>&& _points)
        : points(std::move(_points)), diagonal(points.size(), 0)
    {
        for (auto p : points) {
            assert(p.size() == points.front().size());
        }
    }

    value_t operator()(const index_t i, const index_t j) const
    {
        assert(i < points.size());
        assert(j < points.size());
        return std::sqrt(std::inner_product(
            points[i].begin(), points[i].end(), points[j].begin(), value_t(),
            std::plus<value_t>(),
            [](value_t u, value_t v) { return (u - v) * (u - v); }));
    }

    size_t size() const { return points.size(); }
};


enum file_format {
    LOWER_DISTANCE_MATRIX,
    UPPER_DISTANCE_MATRIX,
    DISTANCE_MATRIX,
    POINT_CLOUD,
    DIPHA,
    SPARSE,
    BINARY
};

static const uint16_t endian_check(0xff00);
static const bool is_big_endian =
    *reinterpret_cast<const uint8_t*>(&endian_check);

template <typename T>
T read(std::istream& input_stream)
{
    T result;
    char* p = reinterpret_cast<char*>(&result);
    if (input_stream.read(p, sizeof(T)).gcount() != sizeof(T))
        return T();
    if (is_big_endian)
        std::reverse(p, p + sizeof(T));
    return result;
}

template <typename T>
bool next_value_stream(T& s, value_t& next_value)
{
    value_t value;
    std::string tmp;
    bool valid = false;

    if (s >> value || !s.eof()) {
        /* Any input stream fails to convert correctly when Inf or NaN is
         * present, this checks if it encounters this kind of values
         * if stof fails to find any valid convertion, it throws an
         * invalid_argument exception*/
        if (s.fail()) {
            s.clear();
            s >> tmp;
            next_value = std::stof(tmp.c_str());
        } else {
            next_value = value;
        }

        /* This check prevent encountering any zero values in distances matrices
         */
        // if (next_value != 0)
        valid = true;
    }

    s.ignore();
    return valid;
}

euclidean_distance_matrix read_point_cloud(std::istream& input_stream)
{
    auto t3 = high_resolution_clock::now();
    std::vector<std::vector<value_t>> points;

    std::string line;
    value_t value;
    while (std::getline(input_stream, line)) {
        std::vector<value_t> point;
        std::istringstream s(line);
        while (s >> value) {
            point.push_back(value);
            s.ignore();
        }
        if (!point.empty())
            points.push_back(point);
        assert(point.size() == points.front().size());
    }

    euclidean_distance_matrix eucl_dist(std::move(points));
    // index_t n = eucl_dist.size();
    // std::cout << "point cloud with " << n << " points in dimension "
    //           << eucl_dist.points.front().size() << std::endl;

    auto t4 = high_resolution_clock::now();
    auto ms_int = duration_cast<milliseconds>(t4 - t3);
    // std::cout << ms_int.count() << "ms\n";

    return eucl_dist;
}

sparse_distance_matrix read_sparse_distance_matrix(std::istream& input_stream)
{
    std::vector<std::vector<index_diameter_t>> neighbors;
    index_t num_edges = 0;

    std::string line;
    while (std::getline(input_stream, line)) {
        std::istringstream s(line);
        size_t i, j;
        value_t value;
        char tmp;
        s >> i;
        s >> tmp;
        s >> j;
        s >> tmp;
        s >> value;
        if (i != j) {
            neighbors.resize(std::max({neighbors.size(), i + 1, j + 1}));
            neighbors[i].push_back({j, value});
            neighbors[j].push_back({i, value});
            ++num_edges;
        }
    }

    for (size_t i = 0; i < neighbors.size(); ++i)
        std::sort(neighbors[i].begin(), neighbors[i].end());

    return sparse_distance_matrix(std::move(neighbors), num_edges);
}

compressed_lower_distance_matrix
read_lower_distance_matrix(std::istream& input_stream)
{
    std::vector<value_t> distances;
    value_t value;
    while (next_value_stream(input_stream, value)) {
        distances.push_back(value);
    }
    std::vector<value_t> diagonal((1 + std::sqrt(1 + 8 * distances.size())) / 2,
            0);

    return compressed_lower_distance_matrix(std::move(distances), std::move(diagonal));
}

compressed_lower_distance_matrix
read_upper_distance_matrix(std::istream& input_stream)
{
    std::vector<value_t> distances;
    value_t value;
    while (next_value_stream(input_stream, value)) {
        distances.push_back(value);
    }
    std::vector<value_t> diagonal((1 + std::sqrt(1 + 8 * distances.size())) / 2,
            0);

    return compressed_lower_distance_matrix(
        compressed_upper_distance_matrix(std::move(distances), std::move(diagonal)));
}

compressed_lower_distance_matrix
read_distance_matrix(std::istream& input_stream)
{
    std::vector<value_t> distances;

    std::string line;
    value_t value;
    for (int i = 0; std::getline(input_stream, line); ++i) {
        std::istringstream s(line);
        for (int j = 0; j < i && next_value_stream(s, value); ++j) {
            distances.push_back(value);
        }
    }
    std::vector<value_t> diagonal((1 + std::sqrt(1 + 8 * distances.size())) / 2,
            0);

    return compressed_lower_distance_matrix(std::move(distances), std::move(diagonal));
}

compressed_lower_distance_matrix read_file(std::istream& input_stream,
                                           const file_format format)
{
    switch (format) {
    case LOWER_DISTANCE_MATRIX:
        return read_lower_distance_matrix(input_stream);
    case UPPER_DISTANCE_MATRIX:
        return read_upper_distance_matrix(input_stream);
    case DISTANCE_MATRIX:
        return read_distance_matrix(input_stream);
    default:
        return read_point_cloud(input_stream);
    }
}


static void BM_ripser(benchmark::State& st, const char* filename, value_t threshold)
{
    index_t dim_max = st.range(0);
    coefficient_t modulus = st.range(1);
    unsigned num_threads = st.range(2);
    file_format format = POINT_CLOUD;
    // value_t threshold = std::numeric_limits<value_t>::max();
    value_t ratio = 1;
    value_t enclosing_radius = std::numeric_limits<value_t>::infinity();
    // const char* filename = "../../data/sphere_3_192_points.dat";

    /* set num_threads to number of concurrent threads supported by the
     * implementation if not set by the user */
    num_threads =
        !num_threads ? std::thread::hardware_concurrency() : num_threads;

    std::ifstream file_stream(filename);
    if (filename && file_stream.fail()) {
        std::cerr << "couldn't open file " << filename << std::endl;
        return;
    }

    std::unique_ptr<sparse_distance_matrix> dist_sparse;
    std::unique_ptr<compressed_lower_distance_matrix> dist_dense;

    if (threshold < std::numeric_limits<value_t>::max()) {
        dist_sparse = std::make_unique<sparse_distance_matrix>(read_point_cloud(filename ? file_stream : std::cin), threshold);
    } else {
        dist_dense = std::make_unique<compressed_lower_distance_matrix>(read_file(filename ? file_stream : std::cin, format));

        value_t min = std::numeric_limits<value_t>::infinity(),
                max = -std::numeric_limits<value_t>::infinity(),
                max_finite = max;
        int num_edges = 0;

        for (auto d : dist_dense->distances) {
            min = std::min(min, d);
            max = std::max(max, d);
            if (d != std::numeric_limits<value_t>::infinity())
                max_finite = std::max(max_finite, d);
            if (d <= threshold)
                ++num_edges;
        }

        if (threshold == std::numeric_limits<value_t>::max()) {
            for (size_t i = 0; i < dist_dense->size(); ++i) {
                value_t r_i = (i==0) ? dist_dense->diagonal[0] : (*dist_dense)(i, 0);
                for (size_t j = 1; j < dist_dense->size(); ++j)
                    r_i = std::max(r_i, (i==j) ? dist_dense->diagonal[i] : (*dist_dense)(i, j));
                enclosing_radius = std::min(enclosing_radius, r_i);
            }

            if (enclosing_radius == std::numeric_limits<value_t>::infinity())
                threshold = max_finite;

            enclosing_radius = std::min(enclosing_radius, max_finite);
        }
    }

    while (st.KeepRunning())
    {
        // auto t1 = high_resolution_clock::now();
        if (threshold < std::numeric_limits<value_t>::max()) {
            sparse_distance_matrix copy_dist = *dist_sparse;

            // sparse_distance_matrix dist(
            //         read_point_cloud(filename ? file_stream : std::cin), threshold);
            auto tmp = ripser<sparse_distance_matrix>(std::move(copy_dist), dim_max, threshold,
                    modulus, num_threads, false);
            tmp.compute_barcodes();
            // tmp.print_barcodes();
        } else {
            compressed_lower_distance_matrix copy_dist = *dist_dense;
            if (threshold == std::numeric_limits<value_t>::max()) {
                auto tmp = ripser<compressed_lower_distance_matrix>(std::move(copy_dist), dim_max,
                        enclosing_radius, modulus, num_threads, false);
                tmp.compute_barcodes();
                // tmp.print_barcodes();
            } else {

                auto tmp = ripser<sparse_distance_matrix>(
                        sparse_distance_matrix(std::move(copy_dist), threshold), dim_max,
                        threshold, modulus, num_threads, false);
                tmp.compute_barcodes();
                // tmp.print_barcodes();
            }
        }
        // auto t2 = high_resolution_clock::now();
        // auto ms_int = duration_cast<milliseconds>(t2 - t1);
        // std::cout << ms_int.count() << "ms\n";
    }
}

// BENCHMARK(BM_ripser)->Args({2, 2, 2})->Iterations(10);
BENCHMARK_CAPTURE(BM_ripser, "sphere", "data/sphere_3_192_points.dat",
                  std::numeric_limits<value_t>::max())
    ->ArgsProduct({{2}, {2}, benchmark::CreateDenseRange(1, MAX_NUMBER_THREADS, 1)})
    ->Iterations(3)
    ->DisplayAggregatesOnly(true)
    ->ReportAggregatesOnly(true);

BENCHMARK_CAPTURE(BM_ripser, "dragon",
                  "data/dragon_vrip.ply.txt_2000_.txt",
                  std::numeric_limits<value_t>::max())
    ->ArgsProduct({{1}, {2}, benchmark::CreateDenseRange(1, MAX_NUMBER_THREADS, 1)})
    ->Iterations(3)
    ->DisplayAggregatesOnly(true)
    ->ReportAggregatesOnly(true);

BENCHMARK_CAPTURE(BM_ripser, "random",
                  "data/random_point_cloud_50_16_.txt",
                  std::numeric_limits<value_t>::max())
    ->ArgsProduct({{7}, {2}, benchmark::CreateDenseRange(1, MAX_NUMBER_THREADS, 1)})
    ->Iterations(3)
    ->DisplayAggregatesOnly(true)
    ->ReportAggregatesOnly(true);

BENCHMARK_CAPTURE(
    BM_ripser, "fractal",
    "data/fractal_9_5_2_linear_edge_list.txt_1866.1116_point_cloud.txt",
    std::numeric_limits<value_t>::max())
    ->ArgsProduct({{2}, {2}, benchmark::CreateDenseRange(1, MAX_NUMBER_THREADS, 1)})
    ->Iterations(3)
    ->DisplayAggregatesOnly(true)
    ->ReportAggregatesOnly(true);

BENCHMARK_CAPTURE(BM_ripser, "o3_1024", "data/o3_1024.txt", 1.8)
    ->ArgsProduct({{3}, {2}, benchmark::CreateDenseRange(1, MAX_NUMBER_THREADS, 1)})
    ->Iterations(3)
    ->DisplayAggregatesOnly(true)
    ->ReportAggregatesOnly(true);

BENCHMARK_CAPTURE(BM_ripser, "o3_4096", "data/o3_4096.txt", 1.4)
    ->ArgsProduct({{3}, {2}, benchmark::CreateDenseRange(1, MAX_NUMBER_THREADS, 1)})
    ->Iterations(2)
    ->DisplayAggregatesOnly(true)
    ->ReportAggregatesOnly(true);

BENCHMARK_CAPTURE(BM_ripser, "torus",
                  "data/clifford_torus_50000.points.txt", 0.15)
    ->ArgsProduct({{2}, {2}, benchmark::CreateDenseRange(1, MAX_NUMBER_THREADS, 1)})
    ->Iterations(1)
    ->DisplayAggregatesOnly(true)
    ->ReportAggregatesOnly(true);

// Run the benchmark
BENCHMARK_MAIN();
