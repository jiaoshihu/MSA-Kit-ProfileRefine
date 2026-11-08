#define MSA_EXPERIMENT_VARIANT 4
#include "internal.hpp"
#include "msakit/core.hpp"

#include <algorithm>
#include <array>
#include <cstdint>
#include <cmath>
#include <limits>
#include <unordered_set>
#include <string_view>
#include <utility>
#include <vector>

namespace msakit {
namespace {
using Score = std::int64_t;

bool is_residue(char value) {
    return value == 'A' || value == 'C' || value == 'G' || value == 'T' ||
           value == 'U' || value == 'N';
}

Score substitution(char center, char query) {
    if (center == 'N' || query == 'N') return 0;
    if (center == query || (center == 'T' && query == 'U') ||
        (center == 'U' && query == 'T')) return 1;
    return -1;
}

struct Projection {
    std::string at_center;
    std::vector<std::string> insertions;
};

Projection align_pair(std::string_view center, std::string_view query) {
    detail::check_score_bounds(center.size(), query.size());
    const auto height = detail::checked_add(center.size(), 1);
    const auto width = detail::checked_add(query.size(), 1);
    const auto cells = detail::checked_multiply(height, width);
    (void)detail::checked_multiply(cells, sizeof(Score));
    if (cells > std::vector<Score>().max_size()) {
        throw detail::Error(detail::kCoreError, "pairwise matrix exceeds container size");
    }

    // Deliberately use the full, easy-to-read DP table in this reference core.
    // It is local to this pair and is released before processing the next query.
    std::vector<Score> scores(cells);
    const auto cell = [&](std::size_t i, std::size_t j) -> Score& {
        return scores[i * width + j];
    };
    for (std::size_t i = 1; i < height; ++i) cell(i, 0) = -static_cast<Score>(i);
    for (std::size_t j = 1; j < width; ++j) cell(0, j) = -static_cast<Score>(j);
    for (std::size_t i = 1; i < height; ++i) {
        for (std::size_t j = 1; j < width; ++j) {
            cell(i, j) = std::max({
                cell(i - 1, j - 1) + substitution(center[i - 1], query[j - 1]),
                cell(i - 1, j) - 1,
                cell(i, j - 1) - 1});
        }
    }

    std::string aligned_center;
    std::string aligned_query;
    const auto longest_trace = detail::checked_add(center.size(), query.size());
    if (longest_trace > aligned_center.max_size()) {
        throw detail::Error(detail::kCoreError, "pairwise traceback exceeds string size");
    }
    aligned_center.reserve(longest_trace);
    aligned_query.reserve(longest_trace);
    std::size_t i = center.size();
    std::size_t j = query.size();
    while (i != 0 || j != 0) {
        // The order of these branches defines the documented tie rule:
        // consume both, consume center only, then consume query only.
        if (i != 0 && j != 0 && cell(i, j) ==
                cell(i - 1, j - 1) + substitution(center[i - 1], query[j - 1])) {
            aligned_center.push_back(center[--i]);
            aligned_query.push_back(query[--j]);
        } else if (i != 0 && cell(i, j) == cell(i - 1, j) - 1) {
            aligned_center.push_back(center[--i]);
            aligned_query.push_back('-');
        } else if (j != 0 && cell(i, j) == cell(i, j - 1) - 1) {
            aligned_center.push_back('-');
            aligned_query.push_back(query[--j]);
        } else {
            throw detail::Error(detail::kCoreError, "invalid pairwise traceback state");
        }
    }
    std::reverse(aligned_center.begin(), aligned_center.end());
    std::reverse(aligned_query.begin(), aligned_query.end());

    Projection projection;
    projection.at_center.assign(center.size(), '-');
    projection.insertions.resize(height);
    std::size_t consumed = 0;
    for (std::size_t column = 0; column < aligned_center.size(); ++column) {
        if (aligned_center[column] == '-') {
            // Slot s follows exactly s consumed center residues. This also
            // accounts for prefix (0) and suffix (center length) insertions.
            projection.insertions[consumed].push_back(aligned_query[column]);
        } else {
            projection.at_center[consumed++] = aligned_query[column];
        }
    }
    return projection;
}
}  // namespace

// Reference star alignment — replaceable core
// Replace this implementation to use another algorithm. Keep input unchanged,
// preserve row-to-input correspondence, and return only aligned strings.

// Experimental adaptive-star center selection.
// Candidate medoid minimizes average distance to other sequences, estimated by
// unique k-mer Jaccard distance plus a small length-difference penalty.
// Sampling prevents quadratic work for large numbers of input sequences.
// For tiny teaching examples keep original first-row behavior.
using KmerSet = std::unordered_set<std::uint64_t>;

KmerSet kmers(std::string_view seq, unsigned k) {
    KmerSet out;
    if (seq.size() < k) return out;
    std::uint64_t word = 0;
    unsigned valid = 0;
    const std::uint64_t mask = (std::uint64_t(1) << (2 * k)) - 1;
    for (char c : seq) {
        if (c == 'U') c = 'T';
        unsigned v = 0;
        if (c == 'A') v = 0;
        else if (c == 'C') v = 1;
        else if (c == 'G') v = 2;
        else if (c == 'T') v = 3;
        else { valid = 0; word = 0; continue; } // skip N-containing windows
        word = ((word << 2) | v) & mask;
        if (++valid >= k) out.insert(word);
    }
    return out;
}

double pair_distance(const KmerSet& a, const KmerSet& b,
                     std::size_t la, std::size_t lb) {
    // Jaccard similarity on sets (NOT alignment identity).
    const auto& small = a.size() <= b.size() ? a : b;
    const auto& large = a.size() <= b.size() ? b : a;
    std::size_t intersection = 0;
    for (auto k : small) if (large.count(k)) ++intersection;
    const std::size_t un = a.size() + b.size() - intersection;
    const double j = un ? double(intersection) / double(un) : 0.0;
    const double len = double(std::max(la, lb) - std::min(la, lb)) /
                       double(std::max(la, lb));
    return (1.0 - j) + 0.2 * len;
}

std::size_t adaptive_kmer_center(const std::vector<SequenceRecord>& input) {
    if (input.size() <= 2) return 0;
    std::size_t max_length = 0;
    for (const auto& rec : input) max_length = std::max(max_length, rec.sequence.size());
    if (max_length < 8) return 0; // preserve teaching-example layouts
    const unsigned k = max_length < 20 ? 3u : 5u;
    const std::size_t n = input.size();
    // Bounded deterministic subset: up to 48 candidates, 96 reference rows.
    // Sampling is evenly spaced in input order; no randomness or external deps.
    auto sample = [n](std::size_t limit) {
        std::vector<std::size_t> ids;
        const std::size_t m = std::min(n, limit);
        ids.reserve(m);
        for (std::size_t i = 0; i < m; ++i)
            ids.push_back(m == 1 ? 0 : i * (n - 1) / (m - 1));
        return ids;
    };
    const auto candidates = sample(48);
    const auto references = sample(96);
    // Cache only selected k-mer sets to bound the prepass memory.
    std::vector<std::size_t> indices = references;
    for (auto i : candidates)
        if (std::find(indices.begin(), indices.end(), i) == indices.end())
            indices.push_back(i);
    std::vector<KmerSet> sets;
    sets.reserve(indices.size());
    for (auto i : indices) sets.push_back(kmers(input[i].sequence, k));
    auto lookup = [&](std::size_t idx) -> const KmerSet& {
        auto it = std::find(indices.begin(), indices.end(), idx);
        return sets[static_cast<std::size_t>(it - indices.begin())];
    };
    std::size_t best = 0;
    double best_score = std::numeric_limits<double>::infinity();
    for (auto c : candidates) {
        double sum = 0.0;
        for (auto r : references) {
            if (c == r) continue;
            sum += pair_distance(lookup(c), lookup(r),
                                 input[c].sequence.size(), input[r].sequence.size());
        }
        if (sum < best_score) { best_score = sum; best = c; }
    }
    return best;
}


// Graph-consistent insertion-slot refinement; experimental research code.
// Only rearranges gaps within pre-existing insertion slots; residue order stays fixed.
// The score sums pairwise interactions between rows (a dense, local correspondence graph).
// For bounded cost, limit rows, insertion widths, and refined slots.
using Shift = std::pair<std::size_t, std::size_t>; // original row, gap shift
using SlotShifts = std::vector<std::vector<Shift>>;

std::uint64_t next_random(std::uint64_t& state) {
    state ^= state << 13; state ^= state >> 7; state ^= state << 17;
    return state;
}

SlotShifts graph_refine_slots(const std::vector<Projection>& projections,
                             std::size_t center_index,
                             const std::vector<std::size_t>& widths,
                             std::size_t number_of_rows) {
    SlotShifts shifts(widths.size());
    if (number_of_rows > 128) return shifts; // explicit computational budget
    std::size_t slots_done = 0;
    for (std::size_t slot = 0; slot < widths.size(); ++slot) {
        const auto width = widths[slot];
        if (width < 2 || width > 32 || slots_done >= 120) continue;
        std::vector<std::size_t> ids;
        std::vector<std::string_view> seqs;
        for (std::size_t r = 0; r < number_of_rows; ++r) {
            if (r == center_index) continue;
            const auto& ins = projections[r].insertions[slot];
            if (!ins.empty()) { ids.push_back(r); seqs.push_back(ins); }
        }
        if (ids.size() < 2 || ids.size() > 48) continue;
        ++slots_done;
        std::vector<std::size_t> pos(ids.size(), 0);
        const auto pair_score = [&](std::size_t a, std::size_t b, std::size_t pa,
                                    std::size_t pb) -> int {
            int score = 0;
            for (std::size_t c = 0; c < width; ++c) {
                const char x = (c >= pa && c-pa < seqs[a].size()) ? seqs[a][c-pa] : '-';
                const char y = (c >= pb && c-pb < seqs[b].size()) ? seqs[b][c-pb] : '-';
                if (x == '-' && y == '-') continue;
                if (x == '-' || y == '-') score -= 1;
                else if (x == y || (x == 'U' && y == 'T') || (x == 'T' && y == 'U')) score += 2;
                else if (x != 'N' && y != 'N') score -= 1;
            }
            return score;
        };
#if MSA_EXPERIMENT_VARIANT == 3
        const auto score_all = [&](const std::vector<std::size_t>& positions) -> long long {
            long long score = 0;
            for (std::size_t a = 0; a < ids.size(); ++a)
                for (std::size_t b = a+1; b < ids.size(); ++b)
                    score += pair_score(a,b,positions[a],positions[b]);
            return score;
        };
#endif
        // V2: coordinate-wise exact local gap-shift search (three sweeps).
        for (int pass = 0; pass < 3; ++pass) {
            bool changed = false;
            for (std::size_t a = 0; a < ids.size(); ++a) {
                int best = 0;
                std::size_t best_pos = pos[a];
                for (std::size_t b = 0; b < ids.size(); ++b)
                    if (a != b) best += pair_score(a,b,pos[a],pos[b]);
                for (std::size_t p = 0; p+seqs[a].size() <= width; ++p) {
                    int candidate = 0;
                    for (std::size_t b = 0; b < ids.size(); ++b)
                        if (a != b) candidate += pair_score(a,b,p,pos[b]);
                    if (candidate > best) { best = candidate; best_pos = p; }
                }
                if (best_pos != pos[a]) { pos[a] = best_pos; changed = true; }
            }
            if (!changed) break;
        }
#if MSA_EXPERIMENT_VARIANT == 3
        // V3: anneal gap shifts; track best-ever state so final score never
        // underperforms V2 on the local graph objective.
        std::vector<std::size_t> best_positions = pos;
        long long current = score_all(pos), best = current;
        std::uint64_t rng = 0x9e3779b97f4a7c15ULL ^ (slot * 0x517cc1b727220a95ULL);
        constexpr int attempts = 250;
        for (int t = 0; t < attempts; ++t) {
            std::size_t a = next_random(rng) % ids.size();
            std::size_t options = width - seqs[a].size() + 1;
            if (options <= 1) continue;
            std::size_t candidate_shift = next_random(rng) % options;
            if (candidate_shift == pos[a]) continue;
            const auto old_shift = pos[a];
            pos[a] = candidate_shift;
            const long long proposed = score_all(pos);
            // Metropolis rule with a cooling temperature and fixed RNG seed.
            const long long delta = proposed - current;
            const double temperature = 0.5 + 6.0 * double(attempts - t) / attempts;
            const double u = double(next_random(rng) >> 11) / 9007199254740992.0;
            const bool accept = delta >= 0 ||
                (double(delta) / temperature > -50.0 &&
                 u < std::exp(double(delta) / temperature));
            if (accept) {
                current = proposed;
                if (current > best) { best = current; best_positions = pos; }
            } else pos[a] = old_shift;
        }
        pos = std::move(best_positions);
#endif
        for (std::size_t a = 0; a < ids.size(); ++a)
            if (pos[a] != 0) shifts[slot].emplace_back(ids[a], pos[a]);
    }
    return shifts;
}

// Profile-profile refinement: rescoring an existing alignment using
// cross-group sum-of-pairs (SP), while preserving every ungapped sequence.
// Each accepted split is a strict global SP improvement under this objective.
namespace {
constexpr int kProfileGapPenalty = -2;
using ProfileCounts = std::array<std::int64_t, 7>; // A,C,G,T(U),N,other,gap
int profile_letter(char ch) {
    switch (ch) {
    case 'A': return 0;
    case 'C': return 1;
    case 'G': return 2;
    case 'T': case 'U': return 3;
    case 'N': return 4;
    case '-': return 6;
    default: return 5;
    }
}
int profile_pair(int a, int b) {
    if (a == 6 && b == 6) return 0;
    if (a == 6 || b == 6) return kProfileGapPenalty;
    if (a == 4 || b == 4) return 0;
    return a == b ? 2 : -1;
}
using GroupColumns = std::vector<ProfileCounts>;
GroupColumns group_columns(const Alignment& alignment,
                           const std::vector<std::size_t>& group) {
    GroupColumns columns;
    const std::size_t width = alignment.rows.front().size();
    columns.reserve(width);
    for (std::size_t c = 0; c < width; ++c) {
        ProfileCounts counts{};
        for (auto r : group) ++counts[profile_letter(alignment.rows[r][c])];
        if (counts[6] != static_cast<std::int64_t>(group.size()))
            columns.push_back(counts); // drop columns empty within group
    }
    return columns;
}
std::vector<std::string> group_alignment(const Alignment& alignment,
                                         const std::vector<std::size_t>& group) {
    std::vector<std::string> result(group.size());
    const std::size_t width = alignment.rows.front().size();
    for (std::size_t c = 0; c < width; ++c) {
        bool any = false;
        for (auto r : group) if (alignment.rows[r][c] != '-') { any = true; break; }
        if (any) for (std::size_t i = 0; i < group.size(); ++i)
            result[i].push_back(alignment.rows[group[i]][c]);
    }
    return result;
}
std::int64_t cross_score(const ProfileCounts& a, const ProfileCounts& b) {
    std::int64_t score = 0;
    for (int i = 0; i < 7; ++i)
        for (int j = 0; j < 7; ++j)
            score += a[i] * b[j] * profile_pair(i, j);
    return score;
}
std::int64_t profile_vs_gaps(const ProfileCounts& col, std::size_t other_rows) {
    return -2LL * (static_cast<std::int64_t>(other_rows) - 0LL) *
           (col[0] + col[1] + col[2] + col[3] + col[4] + col[5]);
}
std::int64_t baseline_cross_score(const Alignment& alignment,
                                  const std::vector<std::size_t>& left,
                                  const std::vector<std::size_t>& right) {
    std::int64_t total = 0;
    const std::size_t width = alignment.rows.front().size();
    for (std::size_t c = 0; c < width; ++c) {
        ProfileCounts a{}, b{};
        for (auto r : left) ++a[profile_letter(alignment.rows[r][c])];
        for (auto r : right) ++b[profile_letter(alignment.rows[r][c])];
        total += cross_score(a,b);
    }
    return total;
}
bool refine_profile_partition(Alignment& alignment,
                              const std::vector<std::size_t>& left,
                              const std::vector<std::size_t>& right) {
    if (left.empty() || right.empty()) return false;
    auto a = group_columns(alignment, left);
    auto b = group_columns(alignment, right);
    if (a.empty() || b.empty()) return false;
    const std::size_t m = a.size(), n = b.size();
    // Safety budget: bound DP memory/time in the evaluation environment.
    if (m > 1200 || n > 1200 || (m+1)*(n+1) > 700000) return false;
    const std::int64_t old_score = baseline_cross_score(alignment, left, right);
    const std::size_t stride = n+1;
    std::vector<std::int64_t> dp((m+1)*stride, 0);
    std::vector<std::uint8_t> trace((m+1)*stride, 0);
    for (std::size_t i=1; i<=m; ++i) {
        dp[i*stride] = dp[(i-1)*stride] + profile_vs_gaps(a[i-1],right.size());
        trace[i*stride] = 1;
    }
    for (std::size_t j=1; j<=n; ++j) {
        dp[j] = dp[j-1] + profile_vs_gaps(b[j-1],left.size());
        trace[j] = 2;
    }
    for (std::size_t i=1; i<=m; ++i) {
        const auto gap_a=profile_vs_gaps(a[i-1],right.size());
        for (std::size_t j=1; j<=n; ++j) {
            std::int64_t best=dp[(i-1)*stride+j-1]+cross_score(a[i-1],b[j-1]);
            std::uint8_t move=0;
            const auto up=dp[(i-1)*stride+j]+gap_a;
            if (up>best) {best=up;move=1;}
            const auto lf=dp[i*stride+j-1]+profile_vs_gaps(b[j-1],left.size());
            if (lf>best) {best=lf;move=2;}
            dp[i*stride+j]=best;
            trace[i*stride+j]=move;
        }
    }
    if (dp[m*stride+n] <= old_score) return false;
    auto pa=group_alignment(alignment,left);
    auto pb=group_alignment(alignment,right);
    std::vector<std::uint8_t> path;
    path.reserve(m+n);
    std::size_t i=m,j=n;
    while (i||j) {
        const auto move=trace[i*stride+j];
        path.push_back(move);
        if(move==0) {--i;--j;}
        else if(move==1) --i;
        else --j;
    }
    std::reverse(path.begin(),path.end());
    Alignment candidate;
    candidate.rows.resize(alignment.rows.size());
    std::size_t ia=0,ib=0;
    for(auto move:path) {
        for(std::size_t r=0;r<left.size();++r)
            candidate.rows[left[r]].push_back(move==2 ? '-' : pa[r][ia]);
        for(std::size_t r=0;r<right.size();++r)
            candidate.rows[right[r]].push_back(move==1 ? '-' : pb[r][ib]);
        if(move!=2) ++ia;
        if(move!=1) ++ib;
    }
    alignment=std::move(candidate);
    return true;
}
void profile_iterative_refine(Alignment& alignment) {
    const std::size_t count=alignment.rows.size();
    if (count<3 || count>64 || alignment.rows.front().size()>1200) return;
    // Two passes over deterministic bipartitions: parity, halves, residue classes.
    // The strictly improving SP acceptance rule ensures deterministic termination.
    for (int iteration=0;iteration<2;++iteration) {
        bool any=false;
        for (int mode=0;mode<3;++mode) {
            std::vector<std::size_t> left,right;
            for(std::size_t r=0;r<count;++r) {
                const bool l = mode==0 ? (r%2==0) :
                               mode==1 ? (r<count/2) : (r%3==0);
                (l?left:right).push_back(r);
            }
            if (refine_profile_partition(alignment,left,right)) any=true;
        }
        if (!any) break;
    }
}
}  // namespace

Alignment align(const std::vector<SequenceRecord>& input) {
    if (input.empty()) {
        throw detail::Error(detail::kCoreError, "the alignment core requires at least one record");
    }
    for (const auto& record : input) {
        if (record.sequence.empty() ||
            !std::all_of(record.sequence.begin(), record.sequence.end(), is_residue)) {
            throw detail::Error(detail::kCoreError,
                                "the alignment core requires nonempty cleaned ACGTUN sequences");
        }
    }
    if (input.size() == 1) return Alignment{{input.front().sequence}};

    const std::size_t center_index = adaptive_kmer_center(input);
    const auto& center = input[center_index].sequence;
    const auto slot_count = detail::checked_add(center.size(), 1);
    std::vector<std::size_t> slot_widths(slot_count, 0);
    // Index by ORIGINAL input row so output order never changes.
    std::vector<Projection> projections(input.size());
    for (std::size_t row = 0; row < input.size(); ++row) {
        if (row == center_index) continue;
        auto projection = align_pair(center, input[row].sequence);
        for (std::size_t slot = 0; slot < slot_count; ++slot) {
            slot_widths[slot] = std::max(slot_widths[slot], projection.insertions[slot].size());
        }
        projections[row] = std::move(projection);
    }

    const auto slot_shifts = graph_refine_slots(projections, center_index, slot_widths, input.size());
    auto output_width = center.size();
    for (const auto count : slot_widths) output_width = detail::checked_add(output_width, count);
    (void)detail::checked_multiply(input.size(), output_width);
    if (output_width > std::string().max_size()) {
        throw detail::Error(detail::kCoreError, "alignment row exceeds string size");
    }

    Alignment result;
    result.rows.resize(input.size());
    for (std::size_t row = 0; row < input.size(); ++row) {
        auto& output = result.rows[row];
        output.reserve(output_width);
        const Projection* projection = row == center_index ? nullptr : &projections[row];
        for (std::size_t slot = 0; slot < slot_count; ++slot) {
            std::size_t used = 0;
            if (projection != nullptr) {
                const auto& insertion = projection->insertions[slot];
                std::size_t shift = 0;
                for (const auto& [original_row, offset] : slot_shifts[slot]) {
                    if (original_row == row) { shift = offset; break; }
                }
                output.append(shift, '-');
                output.append(insertion);
                used = shift + insertion.size();
            }
            // Each insertion string is left-aligned within its slot. A longer
            // query determines the width; shorter queries receive right gaps.
            output.append(slot_widths[slot] - used, '-');
            if (slot != center.size()) {
                output.push_back(projection == nullptr ? center[slot] : projection->at_center[slot]);
            }
        }
    }
    profile_iterative_refine(result);
    return result;
}
}  // namespace msakit
