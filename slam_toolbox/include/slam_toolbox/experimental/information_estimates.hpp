#ifndef SLAM_TOOLBOX__EXPERIMENTAL__THEORETIC_INFORMATION_HPP_
#define SLAM_TOOLBOX__EXPERIMENTAL__THEORETIC_INFORMATION_HPP_

#include "utils.hpp"

class InformationEstimates
{
    typedef std::tuple<int, int, int> map_tuple;

public:
    InformationEstimates(kt_double sensor_range, kt_double resolution, kt_double lambda, kt_double nu);
    InformationEstimates();
    virtual ~InformationEstimates() {}

public:
    // Main function
    std::vector<kt_double> InformationEstimates::findLeastInformativeLaser(std::vector<karto::LocalizedRangeScan*> const& range_scans);
    kt_double findMutualInfo(std::vector<karto::LocalizedRangeScan*> const& range_scans);

private:
    // Mutual information
    kt_double calculateInformationContent(kt_double prob);
    kt_double measurementOutcomeEntropy(map_tuple const& meas_outcome);
    kt_double calculateProbabilityFromLogOdds(kt_double log);
    kt_double mutualInformationFromScans(std::vector<karto::LocalizedRangeScan*> const& range_scans, bool ignore_scan=false, int scan_idx=0);

    // Measurement outcomes probabilities
    void appendCellProbabilities(std::vector<kt_double>& measurements, karto::Vector2<int> const & cell);
    std::unordered_map<map_tuple, kt_double, utils::tuple_hash::HashTuple> computeMeasurementOutcomesHistogram(std::vector<std::vector<kt_double>>& meas_outcm);
    void insertMeasurementOutcome(map_tuple tuple, kt_double probability, std::unordered_map<map_tuple, kt_double, utils::tuple_hash::HashTuple>& map);

    // Measurements calculations <P(free), P(Occ), P(Unk)>
    kt_double calculateScanMassProbabilityBetween(kt_double range_1, kt_double range_2);

private:
    // Data structures
    std::map<karto::Vector2<int>, std::vector<std::vector<kt_double>>> m_cell_probabilities;

    const kt_double l_free = log(0.3 / (1.0 - 0.3));
    const kt_double l_occ = log(0.7 / (1.0 - 0.7));
    const kt_double l_o = log(0.5 / (1.0 - 0.5));

    kt_double m_max_sensor_range;
    kt_double m_cell_resol;
    kt_double m_obs_lambda;
    kt_double m_obs_nu;
    kt_double m_map_dist;
    int m_num_cells;

    // Map grids
    Eigen::MatrixXd m_mutual_grid;
    Eigen::MatrixXi m_visited_grid;

    kt_double m_low_x;
    kt_double m_low_y;

    kt_double m_high_x;
    kt_double m_high_y;
};

#endif
