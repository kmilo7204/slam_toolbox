#include <math.h>
#include <cmath>
#include "slam_toolbox/experimental/information_estimates.hpp"

#include <iostream>

InformationEstimates::InformationEstimates(kt_double sensor_range, kt_double resolution, kt_double lambda, kt_double nu)
: m_max_sensor_range{ sensor_range },
m_cell_resol{ resolution },
m_obs_lambda { lambda },
m_obs_nu { nu }
{
}

InformationEstimates::InformationEstimates()
: m_max_sensor_range{ 25.0 },
m_cell_resol{ 0.1 },
m_obs_lambda{ 0.35 },
m_obs_nu{ 0.28 }
{
}

kt_double InformationEstimates::findMutualInfo(std::vector<karto::LocalizedRangeScan*> const& range_scans)
{
    // double robot_x = 5.0;
    // double robot_y = 5.0;
    // double cell_x = 2.5;
    // double cell_y = 2.0;

    // kt_double angle = atan2(robot_y - cell_y, robot_x - cell_x);

    // std::cout << "Angle: " << angle * (180.0 / M_PI) << std::endl;

    /*
        1. Parate en una celda de la grid inicial.
        2. Identifica cuales laser scans del grupo actual pueden ver esta celda. (Maximum range radius)
        3. Cuales haz de luz del scan seleccionado puede ver esa celda.
        4. Resto de operaciones.
    */

    /*
    Para el punto 0:
        - Debe existir un loop exterior para ir removiendo un laser scan a la vez, con esto identifico cual es el scan
        que aporta menor informacion del grupo dado.
    */

    /*
    Para el punto 1 y 2:
        - Toma el pose del scan actual.
        - Encuentra la distancia del centro del scan con respecto al centro de la celda que nos interesa.
        - Si esa diferencia es menor que el rango del scan, entonces consideramos este scan. (Agregamos ese scan a un vector)

    Esto puedo afrontarlo de forma local, lo que significa que por cada set de scans que reciba,
    puedo encontrar los extremos para poder reescalar el grid, esto potencialmente reduce el tiempo de busqueda.
    */

    /*
    Para el punto 3:
        - Con las lecturas de ese scan, encuentro el haz de luz que este dentro de la celda.
        - Potencialmente puedo encontrar dos lectura, pero por ahora elegire el primero que encuentre.
    */

    // ------------ Start of Grid resizing ------------ //
    // I took the first scan for reference. From this one I extract the min an max points (Poses)
    // std::cout << "-----------------> Step 0 ... " << std::endl;

    m_low_x = range_scans[0]->GetCorrectedPose().GetX();
    m_low_y = range_scans[0]->GetCorrectedPose().GetY();

    m_high_x = range_scans[0]->GetCorrectedPose().GetX();
    m_high_y = range_scans[0]->GetCorrectedPose().GetY();

    // std::cout << "-----------------> Step 1 ... " << std::endl;
    //  Iterating through the group of scans to rescale the grid
    for (const auto & scan : range_scans)
    {
        if (scan == nullptr)
        {
            continue;
        }

        karto::Pose2 pose = scan->GetCorrectedPose();
        karto::PointVectorDouble laser_readings = scan->GetPointReadings(false);
        // std::cout << "Laser reading size: " << laser_readings.size() << std::endl;
        // LaserRangeFinder* pLaserRangeFinder = scan->GetLaserRangeFinder();
        // kt_double rangeThreshold = scan->GetLaserRangeFinder()->GetRangeThreshold();

        kt_double minimumAngle = scan->GetLaserRangeFinder()->GetMinimumAngle();
        kt_double angularResolution = scan->GetLaserRangeFinder()->GetAngularResolution();

        kt_double rangeReading = scan->GetRangeReadings()[0];
        int size = scan->GetRangeReadingsVector().size();
        // std::cout << "Range reading: " << rangeReading << std::endl;
        // std::cout << "Rangescan size: " << size << std::endl;

        // kt_double rangeReading = GetRangeReadings()[i];
        // std::cout << "Sensor pose: " << scan->GetSensorPose().GetHeading() << std::endl;
        // std::cout << "Laser reading 0 - X : " << laser_readings[0].GetX() << std::endl;
        // std::cout << "Laser reading 0 - Y : " << laser_readings[0].GetY() << std::endl;
        // std::cout << "Laser reading 0 - Distance : " << laser_readings[0].Length() << std::endl;
        // // std::cout << "Laser reading 0 - Distance : " << laser_readings[0]. << std::endl;

        // std::cout << "Processing pose: " << pose.GetX() << ", " << pose.GetY() << std::endl;

        // Finding the lower limit pose
        m_low_x = std::min(pose.GetX(), m_low_x);
        m_low_y = std::min(pose.GetY(), m_low_y);

        // Finding the higher limit pose
        m_high_x = std::max(pose.GetX(), m_high_x);
        m_high_y = std::max(pose.GetY(), m_high_y);
    }

    // std::cout << "-----------------> Step 2 ... " << std::endl;


    // std::cout << "Low X Position: " << m_low_x << std::endl;
    // std::cout << "High X Position: " << m_high_x << std::endl;

    // std::cout << "Low Y Position: " << m_low_y << std::endl;
    // std::cout << "High Y Position: " << m_high_y << std::endl;

    // Map dimensions
    kt_double dim_x = std::fabs(m_high_x - m_low_x) + (2 * m_max_sensor_range);
    kt_double dim_y = std::fabs(m_high_y - m_low_y) + (2 * m_max_sensor_range);

    // std::cout << "X Dimension: " << dim_x << std::endl;
    // std::cout << "Y Dimension: " << dim_y << std::endl;

    // Get the number of cells
    int n_cells_x = static_cast<int>(dim_x / m_cell_resol);
    int n_cells_y = static_cast<int>(dim_y / m_cell_resol);

    // Resize the grid with new number of cells
    m_mutual_grid.resize(n_cells_x, n_cells_y);
    m_visited_grid.resize(n_cells_x, n_cells_y);

    // // ------------ End of Grid resizing ------------ //

    // // I would need to move it N times
    // // Cause here I do remove one scan at a time and call the following loop

    // /*
    //     Para encontrar el scan que menor mutual information proporciona debo tener un loop para ir removiendo uno a uno.
    //         1. El loop externo es solamente un for que itera N veces. (Dado por al longitud del vector de scans)
    //         2 Adentro solamente existe un condicional que se salta ese scan en la busqueda actual
    // */
    // // Iterating through the cells
    // kt_double mutual_information = 0.0f;

    // for (int n = 0; n < range_scans.size(); ++n)
    for (int n = 0; n < 1; ++n)
    {
        karto::Pose2 robot_pose_raw = range_scans[n]->GetCorrectedPose();
        // Get the current robot pose and extract the extra range for locating it at the lowest X and Y
        karto::Pose2 local_grid_robot_pose{ robot_pose_raw.GetX() - (m_low_x - m_max_sensor_range), robot_pose_raw.GetY() - (m_low_y - m_max_sensor_range), robot_pose_raw.GetHeading() };
        // std::cout << "Local robot pose: " << local_grid_robot_pose.GetX() << ", " << local_grid_robot_pose.GetY() << std::endl;

        // @NIT
        // This is an ideal point to add an assertion, since the robotr pose must be positive in all of the scenarios

        // Minimum X point
        kt_double lower_limit_x = std::max(0.0, local_grid_robot_pose.GetX() - m_max_sensor_range);
        kt_double upper_limit_x = std::min(dim_x, local_grid_robot_pose.GetX() + m_max_sensor_range);
        kt_double lower_limit_y = std::max(0.0, local_grid_robot_pose.GetY() - m_max_sensor_range);
        kt_double upper_limit_y = std::min(dim_y, local_grid_robot_pose.GetY() + m_max_sensor_range);

        karto::Vector2<kt_double> lower_limit{ lower_limit_x, lower_limit_y };
        karto::Vector2<kt_double> upper_limit{ upper_limit_x, upper_limit_y };

        // @NIT
        // I can use an rvalue reference here
        karto::Vector2<int> lower_limit_cell = utils::grid_operations::getGridPosition(lower_limit, m_cell_resol);
        karto::Vector2<int> upper_limit_cell = utils::grid_operations::getGridPosition(upper_limit, m_cell_resol);

        karto::Vector2<int> local_robot_cell = utils::grid_operations::getGridPosition(local_grid_robot_pose.GetPosition(), m_cell_resol);

        std::cout << " =================== " << std::endl;
        std::cout << "Lower limit cell: " << lower_limit_cell.GetX() << ", " << lower_limit_cell.GetY() << std::endl;
        std::cout << "Upper limit cell: " << upper_limit_cell.GetX() << ", " << upper_limit_cell.GetY() << std::endl;

        // ======================= START OF CELLS ITERATIONS ======================= //
        // for (int i = lower_limit_cell.GetX(); i < local_robot_cell.GetX(); ++i)
        for (int i = lower_limit_cell.GetX(); i < upper_limit_cell.GetX(); ++i)
        {
            for(int j = lower_limit_cell.GetY(); j < upper_limit_cell.GetY(); ++j)
            {
                // This is not fully necessary, but I will leave it here for now ===> Testing
                if ((i < lower_limit_cell.GetX() || i > upper_limit_cell.GetX()) || (j < lower_limit_cell.GetY() || j > upper_limit_cell.GetY()))
                {
                    std::cout << "Out of limits" << std::endl;
                    continue;
                }

                // This one will give me the information
                LaserRangeFinder* laserRangeFinder = range_scans[n]->GetLaserRangeFinder();

                // These will give me the information I need for every laser
                kt_double rangeThreshold = pLaserRangeFinder->GetRangeThreshold();
                kt_double minimumAngle = pLaserRangeFinder->GetMinimumAngle();
                kt_double angularResolution = pLaserRangeFinder->GetAngularResolution();

                // There should be an operator for this, but this is working good.
                karto::Vector2<kt_double> cell_center{ i * m_cell_resol + (m_cell_resol / 2), j * m_cell_resol + (m_cell_resol / 2) };
                kt_double pose_to_cell_angle = atan2(cell_center.GetY() - local_grid_robot_pose.GetY(), cell_center.GetX() - local_grid_robot_pose.GetX());
                pose_to_cell_angle = pose_to_cell_angle > 0 ? pose_to_cell_angle : (2.0 * M_PI + pose_to_cell_angle);

                kt_double robot_heading = robot_pose_raw.GetHeading();
                robot_heading = robot_heading < 0.0 ? (2.0 * M_PI + robot_heading) : robot_heading;

                int laser_index = robot_heading / 0.0174533; // This one will be changed

                // So remember that now I need to count from the minimum range.
                // There should be an operator for this, but this is working good.
                karto::Vector2<kt_double> cell_center{ i * m_cell_resol + (m_cell_resol / 2), j * m_cell_resol + (m_cell_resol / 2) };
                kt_double pose_to_cell_angle = atan2(cell_center.GetY() - local_grid_robot_pose.GetY(), cell_center.GetX() - local_grid_robot_pose.GetX());
                pose_to_cell_angle = pose_to_cell_angle > 0 ? pose_to_cell_angle : (2.0 * M_PI + pose_to_cell_angle);

                kt_double robot_heading = robot_pose_raw.GetHeading();
                robot_heading = robot_heading < 0.0 ? (2.0 * M_PI + robot_heading) : robot_heading;

                int laser_index = robot_heading / 0.0174533;
                int angle_idx = pose_to_cell_angle / 0.0174533;
                int aux_index = angle_idx;
                if (laser_index > 180)
                {
                    aux_index = 360 - laser_index + angle_idx;
                }
                else if (laser_index >= 0 && laser_index <= 180)
                {
                    aux_index = laser_index + (laser_index - angle_idx);
                }

                kt_double reading_distance_raw = range_scans[n]->GetRangeReadings()[aux_index];
                kt_double reading_distance = range_scans[n]->GetRangeReadings()[aux_index];

                if (reading_distance > m_max_sensor_range)
                {
                    continue;
                }

                kt_double distance_to_cell = sqrt(pow(local_grid_robot_pose.GetX() - cell_center.GetX(), 2) + pow(local_grid_robot_pose.GetY() - cell_center.GetY(), 2));
                if (reading_distance < distance_to_cell)
                {
                    continue;
                }
                else
                {
                    reading_distance = distance_to_cell;
                }


                kt_double point_x = local_grid_robot_pose.GetX() + (reading_distance * cos(pose_to_cell_angle));
                kt_double point_y = local_grid_robot_pose.GetY() + (reading_distance * sin(pose_to_cell_angle));

                karto::Vector2<int> laser_point = utils::grid_operations::getGridPosition({point_x, point_y}, m_cell_resol);

                if( sqrt(pow((i * m_cell_resol) - point_x, 2) + pow((j * m_cell_resol - point_y), 2)) > 0.5)
                {
                    // Aqui probablmente pueda correr un for con los siguientes con dos indices anteriores y dos indices siguientes
                    continue;
                }

                std::cout << "Robot heading: " << robot_heading * (180 / M_PI) << std::endl;
                std::cout << "Angle to cell: " << pose_to_cell_angle * (180 / M_PI) << std::endl;
                std::cout << "Laser index: " << angle_idx << std::endl;
                std::cout << "Guessed index: " << aux_index << std::endl;
                std::cout << "Guessed angle: " << aux_index * (180 / M_PI) << std::endl;
                std::cout << "Distance to cell: " << distance_to_cell << std::endl;
                std::cout << "Reading distance: " << reading_distance << std::endl;
                std::cout << "Reading distance raw: " << reading_distance_raw << std::endl;
                std::cout << "Robot pose: " << local_robot_cell.GetX() << ", " << local_robot_cell.GetY() << std::endl;
                std::cout << "Laser coordinates: " << point_x << ", " << point_y << std::endl;
                std::cout << "Cell center: " << cell_center.GetX() << ", " << cell_center.GetY() << std::endl;
                std::cout << "Laser cell: " << laser_point.GetX() << ", " << laser_point.GetY() << std::endl;
                std::cout << "Current cell: " << i << ", " << j << std::endl;
                std::cout << " =================== " << std::endl;

                // Inidividual cell limits
                kt_double limit_x = i * m_cell_resol;
                kt_double limit_y = j * m_cell_resol;

                std::pair<std::vector<kt_double>, std::vector<kt_double>> intersections = utils::grid_operations::computeLineBoxIntersection(
                robot_pose.GetPosition(), transformed_laser, robot_grid_position, beam_grid, limit_x, limit_y, m_cell_resol);

                // In case we did not find a measurement
                if (intersections.first.size() == 0)
                {
                    continue;
                }

                // Enter (d1) and Exit (d2) distances
                std::vector<kt_double> distances;
                for (int k = 0; k < intersections.first.size(); ++k)
                {
                    // From robot position to intersection points
                    karto::Vector2<kt_double> intersection{intersections.first[k], intersections.second[k]};
                    kt_double distance = robot_pose.GetPosition().Distance(intersection);
                    distances.push_back(distance);
                }

                // Measurement outcomes vector {Pfree, Pocc, Pun}
                std::vector<kt_double> probabilities {
                    calculateScanMassProbabilityBetween(distances[1], m_max_sensor_range),
                    calculateScanMassProbabilityBetween(distances[0], distances[1]),
                    calculateScanMassProbabilityBetween(0.0f, distances[0])
                };

            }
        }
        std::cout << "<------------------------->" << std::endl;
    }
    return 5.0;
}

std::vector<kt_double> InformationEstimates::findLeastInformativeLaser(std::vector<karto::LocalizedRangeScan*> const& range_scans)
{
    /**
     * Find and return the Laser Scan index which provide the less mutual information in a set of Laser Scans
     * Arguments:
        * range_scans [std::vector<karto::LocalizedRangeScan*>]: Vector of LocalizedRangeScan pointers
     * Return:
        * std::tuple<int, kt_double>: Tuple containing the index of the LozalizedRangeScan and its corresponding mutual information
    */
    // std::cout << "Longitude: " << range_scans.size() << std::endl;

    // This is under testing
    m_low_x = range_scans[0]->GetCorrectedPose().GetX();
    m_low_y = range_scans[0]->GetCorrectedPose().GetY();

    m_high_x = range_scans[0]->GetCorrectedPose().GetX();
    m_high_y = range_scans[0]->GetCorrectedPose().GetY();

    for (const auto & scan : range_scans)
    {
        if (scan == nullptr)
        {
            continue;
        }

        karto::Pose2 pose = scan->GetCorrectedPose();
        // std::cout << "Processing pose: " << pose.GetX() << ", " << pose.GetY() << std::endl;

        // Finding the closest pose
        m_low_x = pose.GetX() < m_low_x ? pose.GetX() : m_low_x;
        m_low_y = pose.GetY() < m_low_y ? pose.GetY() : m_low_y;

        // Finding the farthest pose
        m_high_x = pose.GetX() > m_high_x ? pose.GetX() : m_high_x;
        m_high_y = pose.GetY() > m_high_y ? pose.GetY() : m_high_y;
    }

    // Margins for the closest points
    m_low_x -= m_max_sensor_range;
    m_low_y -= m_max_sensor_range;

    // Margins for the farthest points
    m_high_x += m_max_sensor_range;
    m_high_y += m_max_sensor_range;

    // ===== Grid resizing =====
    // Map dimensions
    kt_double dist_x = std::fabs(m_high_x) + std::fabs(m_low_x);
    kt_double dist_y = std::fabs(m_high_y) + std::fabs(m_low_y);

    // std::cout << "X: " << m_low_x << ", " << m_high_x << std::endl;
    // std::cout << "Y: " << m_low_y << ", " << m_high_y << std::endl;
    // std::cout << "Distances: " << dist_x << ", " << dist_y << std::endl;

    // Get the number of cells
    int n_cells_x = static_cast<int>(dist_x / m_cell_resol);
    int n_cells_y = static_cast<int>(dist_y / m_cell_resol);
    // std::cout << "Number of cells: " << n_cells_x << ", " << n_cells_y << std::endl;

    // Resize the grid with new number of cells
    m_mutual_grid.resize(n_cells_x, n_cells_y);
    m_visited_grid.resize(n_cells_x, n_cells_y);
    // // ===== Grid resizing =====

    std::vector<kt_double> scans_mut_info;
    scans_mut_info.reserve(range_scans.size());


    // Total mutual information
    kt_double map_mut_info = mutualInformationFromScans(range_scans);
    // Clearing the cells for the next time it is called
    m_mutual_grid.setZero();

    for (int i = 0; i < range_scans.size(); ++i)
    {
        // Find mutual information when extracting one laser scan
        kt_double temp_mut_info = mutualInformationFromScans(range_scans, true, i);
        scans_mut_info.emplace_back(map_mut_info - temp_mut_info);
        // Clearing the cells for the next time it is called
        m_mutual_grid.setZero();
    }

    // Finding the less informative laser scan
    std::vector<kt_double>::iterator it_min;
    it_min = std::min_element(scans_mut_info.begin(), scans_mut_info.end());
    int idx = it_min - scans_mut_info.begin();

    return scans_mut_info;
}

kt_double InformationEstimates::mutualInformationFromScans(std::vector<karto::LocalizedRangeScan*> const& range_scans, bool ignore_scan, int scan_idx)
{
    /**
     * Append measured porbabilities for the given cell
     * Arguments:
        * range_scans [std::vector<karto::LocalizedRangeScan*>]: Vector of LocalizedRangeScan pointers
        * ignore_scan [bool]: Indicates the type of calculation
        * scan_idx [int]: Index within the set to be ignored
     * Return:
        * kt_double
    */

    int curr_idx = 0;
    std::vector<karto::Vector2<int>> vis_cells;

    for (auto & scan : range_scans)
    {
        if ((curr_idx == scan_idx && ignore_scan) || (scan == nullptr))
        {
            ++curr_idx;
            continue;
        }

        karto::Pose2 robot_pose_raw = scan->GetCorrectedPose();
        karto::Pose2 robot_pose(robot_pose_raw.GetX() - m_low_x, robot_pose_raw.GetY() - m_low_y, robot_pose_raw.GetHeading());
        karto::PointVectorDouble laser_readings = scan->GetPointReadings(true);
        karto::Vector2<int> robot_grid = utils::grid_operations::getGridPosition(robot_pose.GetPosition(), m_cell_resol);

        // Set as false the current boolean map
        m_visited_grid.setZero();

        for (int i = 0; i < laser_readings.size(); ++i)
        {
            karto::Vector2<kt_double> transformed_laser{laser_readings[i].GetX() - m_low_x, laser_readings[i].GetY() - m_low_y};

            // Laser final cell
            karto::Vector2<int> beam_grid = utils::grid_operations::getGridPosition(transformed_laser, m_cell_resol);

            // Visited cells by this laser beam
            std::vector<karto::Vector2<int>> cells = utils::grid_operations::rayCasting(robot_grid, beam_grid);

            for (auto & cell : cells)
            {
                vis_cells.push_back(cell);

                // Inidividual cell limits
                kt_double limit_x = cell.GetX() * m_cell_resol;
                kt_double limit_y = cell.GetY() * m_cell_resol;

                std::pair<std::vector<kt_double>, std::vector<kt_double>> intersections = utils::grid_operations::computeLineBoxIntersection(
                    robot_pose.GetPosition(), transformed_laser, robot_grid, beam_grid, limit_x, limit_y, m_cell_resol);

                if (intersections.first.size() == 0)
                {
                    continue;
                }

                // Enter (d1) and Exit (d2) distances
                std::vector<kt_double> distances;
                for (int k = 0; k < intersections.first.size(); ++k)
                {
                    // From robot position to intersection points
                    karto::Vector2<kt_double> intersection{intersections.first[k], intersections.second[k]};
                    kt_double distance = robot_pose.GetPosition().Distance(intersection);
                    distances.push_back(distance);
                }

                // Measurement outcomes vector {Pfree, Pocc, Pun}
                std::vector<kt_double> probabilities {
                    calculateScanMassProbabilityBetween(distances[1], m_max_sensor_range),
                    calculateScanMassProbabilityBetween(distances[0], distances[1]),
                    calculateScanMassProbabilityBetween(0.0f, distances[0])
                };

                // Appending new measurement outcomes for the current cell
                appendCellProbabilities(probabilities, cell);
            }
        }
        ++curr_idx;
    }

    /*
        - Changing the calculation point to here have one benefit.
            - I will only calculate the mutual information once.

        - Calculate the Measurement outcomes histograms for each cell.

        - I need to track which cells where actually observed from one of the lasers.
    */
    for (auto & cell : vis_cells)
    {
        std::unordered_map<map_tuple, kt_double, utils::tuple_hash::HashTuple> meas_out_prob = computeMeasurementOutcomesHistogram(m_cell_probabilities.at(cell));
        kt_double cell_mutual_inf = 0.0f;
        for (auto& pair : meas_out_prob)
        {
            cell_mutual_inf +=  pair.second * measurementOutcomeEntropy(pair.first);
        }

        // Mutual information of cell x, y given a set of measurements
        m_mutual_grid(cell.GetX(), cell.GetY()) = 1.0 - cell_mutual_inf;
    }

    m_cell_probabilities.clear();

    return m_mutual_grid.sum();
}


void InformationEstimates::appendCellProbabilities(std::vector<kt_double>& measurements, karto::Vector2<int> const & cell)
{
    /**
     * Append measured porbabilities for the given cell
     * Arguments:
        * measurements [std::vector<kt_double>]: Vector of LocalizedRangeScan pointers
        * cell [karto::Vector2<int>]: Cell for appending the data
     * Return:
        * Void
    */

    /*
        This is what is going on:
            - In mutualInformationFromScans is a line for cleaning the visited cells m_visited_grid.setZero(), it takes place
            right before we start processing the visited cells for a new scan.
            - Some operations are done to calculate the probability of see a cell as occupied, free or unknown.
            - The probabilities are given to appendCellProbabilities, so we store them into our map.
            - Inside appendCellProbabilities we evaluate if the each cell in the current laser reading has been visited:
                - Since we consider that a cell can be seen only once for one beam of the laser scan. In theory two beams should not
                see the same cell, but as this is a grid map, it might happen due to the resolution. For avoiding it, we are performing
                the following evaluation on each of the laser beams of the current laser scan.
                    - Has this cell been seen by a beam?. If not then add it to our map, and mark this cells as a visited. EOF the function
                    - Is the cell in the map and has been marked as visited by a beam?. If yes, then replace it because we don't want to have
                    duplicate values
                    - Is the cell in the map and has not been marked as visited by the current beam.
                    Case 3 Other beam of other laser scan sees the cell and it has not been marked as visited, because this new laser
                    just discovered that cell

        I need to triger some testing for this case
    */
    std::map<karto::Vector2<int>, std::vector<std::vector<kt_double>>>::iterator it_cell;
    it_cell = m_cell_probabilities.find(cell);

    /*
        Llego la medicion de una celda especifica

        1. Si esa celda no esta presente, entonces agregamos la celda al mapa de probabilidades
        2. La celda esta presente en el mapa de probabilidades y fue marcada como visitada,
            entonces solo reemplazo las probabilidades.
        3. La celda no ha sido marcada como visitada,
    */
    if (it_cell == m_cell_probabilities.end())
    {
        // std::cout << "Case 1: " << cell.GetX() << ", " << cell.GetY() << std::endl;
        // Cell is not present in the map, so append it
        m_cell_probabilities.insert(std::pair<karto::Vector2<int>, std::vector<std::vector<kt_double>>>(
            cell, {measurements}));
        m_visited_grid(cell.GetX(), cell.GetY()) = 1;
    }
    else
    {
        if(m_visited_grid(cell.GetX(), cell.GetY()) == 1) // It must be the longitude
        {
            // Compare the unknown probability, the smallest it is the most information we will have
            // from the occupied or free state
            // std::cout << "Case 2: " << cell.GetX() << ", " << cell.GetY() << std::endl;
            int idx = it_cell->second.size() - 1;
            if(measurements[2] < it_cell->second[idx][2])
            {
                // Replacing
                it_cell->second[idx][0] = measurements[0];
                it_cell->second[idx][1] = measurements[1];
                it_cell->second[idx][2] = measurements[2];
            }
        }
        else
        {
            /*
                Si la celda ya esta agregada, entonces debo hacer un append de las probabilidades
            */
            // Cell is already in the map, only add the next measurement outcome
            // it_cell->second.push_back({measurements[0], measurements[1], measurements[2]});

            // std::cout << "Visiting this code here" << std::endl;
            // std::cout << "Case 3: " << cell.GetX() << ", " << cell.GetY() << std::endl;
            it_cell->second.push_back(measurements);
            m_visited_grid(cell.GetX(), cell.GetY()) = 1;
        }
    }
}

kt_double InformationEstimates::calculateInformationContent(kt_double prob)
{
    /**
     * Calculate the information content or self-information based on the probability of cell being occupied
     * Arguments:
        * prob [kt_double]: Probability of being occupied
     * Return:
        * kt_double: Information content
    */
    return - (prob * log2(prob)) -  ((1 - prob) * log2(1 - prob));
}

kt_double InformationEstimates::calculateProbabilityFromLogOdds(kt_double log)
{
    /**
     * Map Log-Odds into probability
     * Arguments:
        * log [kt_double]: Log-Odds
     * Return:
        * kt_double: Probability
    */
    return (exp(log) / (1 + exp(log)));
}

void InformationEstimates::insertMeasurementOutcome(map_tuple tuple, kt_double probability, std::unordered_map<map_tuple, kt_double, utils::tuple_hash::HashTuple>& map)
{
    int key_free, key_occ, key_unk;
    std::tie(key_free, key_occ, key_unk) = tuple;

    if (map[std::make_tuple(key_free, key_occ, key_unk)])
    {
        map[std::make_tuple(key_free, key_occ, key_unk)] += probability;
    }
    else
    {
        map[std::make_tuple(key_free, key_occ, key_unk)] = probability;
    }
}

std::unordered_map<InformationEstimates::map_tuple, kt_double, utils::tuple_hash::HashTuple> InformationEstimates::computeMeasurementOutcomesHistogram(std::vector<std::vector<kt_double>>& meas_outcm)
{
    /**
     * Implementation of algorithm 1
     * Compute all the possible combinations of a grid cell, given a set of measurement outcomes
     * Arguments:
        * meas_outcm [std::vector<std::vector<kt_double>>]: Vector of measurement outcomes in the form {p_free, p_occ, p_unk}
     * Return:
        * std::unordered_map<InformationEstimates::map_tuple, kt_double, utils::tuple_hash::HashTuple>: Map of combination, it contains the combination and its probability
    */

    std::unordered_map<map_tuple, kt_double, utils::tuple_hash::HashTuple> probabilities_map;
    std::unordered_map<map_tuple, kt_double, utils::tuple_hash::HashTuple> temp_map;

    // Clear the temporal map
    temp_map.clear();
    probabilities_map.clear();

    // Root
    if (meas_outcm.size() == 0)
    {
        probabilities_map[std::make_tuple(0, 0, 0)] = 1.0f;
    }

    // First measurement
    probabilities_map[std::make_tuple(1, 0, 0)] = meas_outcm[0][0]; // Probability free
    probabilities_map[std::make_tuple(0, 1, 0)] = meas_outcm[0][1]; // Probability occupied
    probabilities_map[std::make_tuple(0, 0, 1)] = meas_outcm[0][2]; // Probability unknown

    for (int r = 1; r < meas_outcm.size(); ++r)
    {
        // Measurement outcome probability
        kt_double free_prop = meas_outcm[r][0];
        kt_double occ_prop = meas_outcm[r][1];
        kt_double unk_prop = meas_outcm[r][2];

        // Temporal map to only take the last level combination
        temp_map = probabilities_map;
        probabilities_map.clear();

        for (auto & combination : temp_map)
        {
            int key_free, key_occ, key_unk;
            std::tie(key_free, key_occ, key_unk) = combination.first;

            // Adding next level measurement outcomes
            insertMeasurementOutcome(std::make_tuple(key_free + 1, key_occ, key_unk), combination.second * free_prop, probabilities_map);
            insertMeasurementOutcome(std::make_tuple(key_free, key_occ + 1, key_unk), combination.second * occ_prop, probabilities_map);
            insertMeasurementOutcome(std::make_tuple(key_free, key_occ, key_unk + 1), combination.second * unk_prop, probabilities_map);
        }
        temp_map.clear();
    }
    return probabilities_map;
}

kt_double InformationEstimates::measurementOutcomeEntropy(map_tuple const& meas_outcome)
{
    /**
     * Calculate the measurement outcome entropy
        * Calculate Log-Odds from initial probability guess
        * Calculate the probability from those logs
        * Calculate the entropy with the retrieved probability
     * Arguments:
        * meas_outcome [map_tuple]: Measurement outcome in the form {p_free, p_occ, p_unk}
     * Return:
        * kt_double: Measurement outcome entropy
    */
    int num_free, num_occ, num_unk;
    std::tie(num_free, num_occ, num_unk) = meas_outcome;
    kt_double log_occ = (num_free * l_free) + (num_occ * l_occ) - ((num_free + num_occ - 1) * l_o);
    kt_double prob_occ = calculateProbabilityFromLogOdds(log_occ);
    return calculateInformationContent(prob_occ);
}

kt_double InformationEstimates::calculateScanMassProbabilityBetween(kt_double range_1, kt_double range_2)
{
    /**
     * Calculate the mass probability of a cell being observed by a given measurement
     * Arguments:
        * range_1 [kt_double]: Lower bound
        * range_2 [kt_double]: Upper bound
     * Return:
        * kt_double: Mass probability
    */
    range_1 = (range_1 > m_max_sensor_range) ? m_max_sensor_range : range_1;
    range_2 = (range_2 > m_max_sensor_range) ? m_max_sensor_range : range_2;

    return m_obs_nu * (exp(-m_obs_lambda*range_1) - exp(-m_obs_lambda*range_2));
}

