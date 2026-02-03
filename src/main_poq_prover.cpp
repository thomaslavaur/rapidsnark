#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <stdexcept>
#include <cstdint>
#include <nlohmann/json.hpp>
#include "prover.h"
#include "fileloader.hpp"

using json = nlohmann::json;

static std::vector<u_int32_t> LoadMutableIndices(const std::string &filename)
{
    std::ifstream indicesFile(filename);
    if (!indicesFile) {
        throw std::runtime_error("Unable to open mutable indices file");
    }

    json indicesJson;
    indicesFile >> indicesJson;

    if (!indicesJson.is_array()) {
        throw std::runtime_error("Mutable indices file must contain a JSON array");
    }

    std::vector<u_int32_t> indices;
    for (const auto &item : indicesJson) {
        if (!item.is_number_unsigned()) {
            throw std::runtime_error("Mutable indices must be unsigned integers");
        }
        indices.push_back(item.get<u_int32_t>());
    }

    return indices;
}

int main(int argc, char **argv)
{
    if (argc < 6) {
        std::cerr << "Invalid number of parameters" << std::endl;
        std::cerr << "Usage: poq_prover <circuit.zkey> <base.wtns> <mutable_indices.json> <count>"
                  << " <witness_1.wtns> ... <witness_n.wtns>"
                  << " <proof_1.json> ... <proof_n.json>"
                  << " <public_1.json> ... <public_n.json>" << std::endl;
        return EXIT_FAILURE;
    }

    try {
        const std::string zkeyFilename = argv[1];
        const std::string baseWtnsFilename = argv[2];
        const std::string indicesFilename = argv[3];
        const int count = std::stoi(argv[4]);

        if (count <= 0) {
            throw std::runtime_error("Proof count must be positive");
        }

        const int expectedArgs = 5 + count * 3;
        if (argc != expectedArgs) {
            throw std::runtime_error("Invalid number of parameters for the given proof count");
        }

        std::vector<std::string> witnessFiles;
        std::vector<std::string> proofFiles;
        std::vector<std::string> publicFiles;

        witnessFiles.reserve(count);
        proofFiles.reserve(count);
        publicFiles.reserve(count);

        int offset = 5;
        for (int i = 0; i < count; ++i) {
            witnessFiles.emplace_back(argv[offset + i]);
        }
        offset += count;
        for (int i = 0; i < count; ++i) {
            proofFiles.emplace_back(argv[offset + i]);
        }
        offset += count;
        for (int i = 0; i < count; ++i) {
            publicFiles.emplace_back(argv[offset + i]);
        }

        const auto mutableIndices = LoadMutableIndices(indicesFilename);

        BinFileUtils::FileLoader zkeyFile(zkeyFilename);
        BinFileUtils::FileLoader baseWtnsFile(baseWtnsFilename);

        void *poqProver = nullptr;
        char errorMsg[1024];

        int error = groth16_poq_prover_create(
            &poqProver,
            zkeyFile.dataBuffer(),
            zkeyFile.dataSize(),
            baseWtnsFile.dataBuffer(),
            baseWtnsFile.dataSize(),
            mutableIndices.data(),
            mutableIndices.size(),
            errorMsg,
            sizeof(errorMsg));

        if (error != PROVER_OK) {
            throw std::runtime_error(errorMsg);
        }

        unsigned long long publicSize = 0;
        error = groth16_public_size_for_zkey_buf(
            zkeyFile.dataBuffer(),
            zkeyFile.dataSize(),
            &publicSize,
            errorMsg,
            sizeof(errorMsg));

        if (error != PROVER_OK) {
            groth16_poq_prover_destroy(poqProver);
            throw std::runtime_error(errorMsg);
        }

        unsigned long long proofSize = 0;
        groth16_proof_size(&proofSize);

        std::vector<char> publicBuffer(publicSize);
        std::vector<char> proofBuffer(proofSize);

        for (int i = 0; i < count; ++i) {
            BinFileUtils::FileLoader wtnsFile(witnessFiles[i]);
            unsigned long long currentProofSize = proofSize;
            unsigned long long currentPublicSize = publicSize;

            error = groth16_poq_prover_prove(
                poqProver,
                wtnsFile.dataBuffer(),
                wtnsFile.dataSize(),
                proofBuffer.data(),
                &currentProofSize,
                publicBuffer.data(),
                &currentPublicSize,
                errorMsg,
                sizeof(errorMsg));

            if (error != PROVER_OK) {
                groth16_poq_prover_destroy(poqProver);
                throw std::runtime_error(errorMsg);
            }

            std::ofstream proofFile(proofFiles[i]);
            proofFile.write(proofBuffer.data(), currentProofSize);

            std::ofstream publicFile(publicFiles[i]);
            publicFile.write(publicBuffer.data(), currentPublicSize);
        }

        groth16_poq_prover_destroy(poqProver);

    } catch (std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return EXIT_FAILURE;
    }

    exit(EXIT_SUCCESS);
}
