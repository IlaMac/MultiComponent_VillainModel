#pragma once
#include <mpi.h>



namespace mpi {

    std::vector<std::vector<double>>
    gather(const std::vector<double> &send_data, int source, int destination, int num_ranks, MPI_Datatype mpi_type) {
        std::vector<double> receive_data;
        if (source == destination) {
            receive_data.resize(num_ranks * send_data.size());
        }

        MPI_Gather(send_data.data(), send_data.size(), mpi_type, receive_data.data(),send_data.size(), mpi_type, destination, MPI_COMM_WORLD);
        std::vector<std::vector<double>> split_data;
        if (source == destination) {
            split_data.resize(num_ranks);
            int offset = 0;
            for (auto &vec : split_data) {
                vec.resize(send_data.size());
//                auto last = temp_potential.begin() + offset + vil.potential.size();
//                auto first = temp_potential.begin() + offset ;
//                if(last > temp_potential.size()) throw std::runtime_error(...);
                vec = std::vector<double>(receive_data.begin() + offset,
                                          receive_data.begin() + offset + send_data.size());
                offset += send_data.size();
            }
        }
        return split_data;
    }

//    void scatter(const std::vector<double> &merged_data,std::vector<double> &receive_data, int source, int num_ranks, MPI_Datatype mpi_type) {
//        // Send data comes from root. It is a big vector with all the merged data that is going to be split
//
//        size_t receive_size = merged_data.size() / num_ranks;
//        if(receive_size != receive_data.size())
//            throw std::runtime_error("Incorrect size on receive data. Expected "
//                                     + std::to_string(receive_size) + " got "
//                                     + std::to_string(receive_data.size()));
//
//        MPI_Scatter(merged_data.data(), receive_size, MPI_DOUBLE, receive_data.data(),receive_size , MPI_DOUBLE, source, MPI_COMM_WORLD);
//    }

    void scatter(const std::vector<std::vector<double>> &split_data,std::vector<double> &receive_data, int source, int destination, int num_ranks, MPI_Datatype mpi_type) {
        // Send data comes from root. It is a big vector with all the merged data that is going to be split
        // Start with merging the vecotr into one big vector for scattering
        std::vector<double> merged_data;

        if (source == destination) {
            // Flatten split_data
            merged_data.reserve(receive_data.size()*num_ranks);
            for(auto & s : split_data){
                std::copy(s.begin(),s.end(), std::back_inserter(merged_data));
            }

            size_t receive_size = merged_data.size() / num_ranks;
            if(receive_size != receive_data.size())
                throw std::runtime_error("Incorrect size on receive data. Expected "
                                         + std::to_string(receive_size) + " got "
                                         + std::to_string(receive_data.size()));
        }


        MPI_Scatter(merged_data.data(), receive_data.size(), MPI_DOUBLE, receive_data.data(), receive_data.size(), MPI_DOUBLE, source, MPI_COMM_WORLD);
    }

}




//        if constexpr( is_std_vector(T)){
//
//        }else ( is_pointer(T)){
//                T receive_data(data.size());
//                MPI_Gather(data.data(), 1 , vill_struct, receive_data.data(); 1, vill_struct, PTp.root, MPI_COMM_WORLD);
//                potential = mpi::gather(vil.potential, Ptp.root);
//
//        }

//    }
//
//    template<typename T>
//    T gather(T & data, int root, MPI_Datatype mpi_type){
//        T receive_data(number of temps * data.size());
//        MPI_Gather(data.data(), 1 , vill_struct, receive_data.data(); 1, vill_struct, PTp.root, MPI_COMM_WORLD);
//        potential = mpi::gather(vil.potential, Ptp.root);
//        if constexpr( is_std_vector(T)){
//
//        }else ( is_pointer(T)){
//                T receive_data(data.size());
//                MPI_Gather(data.data(), 1 , vill_struct, receive_data.data(); 1, vill_struct, PTp.root, MPI_COMM_WORLD);
//                potential = mpi::gather(vil.potential, Ptp.root);
//
//        }

//    }

//
//    template<typename T>
//    T scatter(T & data, int root, MPI_Datatype mpi_type){
//        T receive_data(number of temps * data.size());
//        MPI_Gather(data.data(), 1 , vill_struct, receive_data.data(); 1, vill_struct, PTp.root, MPI_COMM_WORLD);
//        potential = mpi::gather(vil.potential, Ptp.root);
//        if constexpr( is_std_vector(T)){
//
//        }else ( is_pointer(T)){
//                T receive_data(data.size());
//                MPI_Gather(data.data(), 1 , vill_struct, receive_data.data(); 1, vill_struct, PTp.root, MPI_COMM_WORLD);
//                potential = mpi::gather(vil.potential, Ptp.root);
//
//        }

//    }

//}