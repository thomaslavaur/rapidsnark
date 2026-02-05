#include <algorithm>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <stdexcept>

#include <alt_bn128.hpp>
#include <nlohmann/json.hpp>

#include "binfile_utils.hpp"
#include "fileloader.hpp"
#include "groth16.hpp"
#include "wtns_utils.hpp"
#include "zkey_utils.hpp"

using json = nlohmann::json;

namespace {

bool PrimeIsValid(mpz_srcptr prime)
{
    mpz_t altBbn128r;

    mpz_init(altBbn128r);
    mpz_set_str(altBbn128r,
                "21888242871839275222246405745257275088548364400416034343698204186575808495617",
                10);

    const bool is_valid = (mpz_cmp(prime, altBbn128r) == 0);

    mpz_clear(altBbn128r);

    return is_valid;
}

std::vector<u_int32_t> BuildModifiedIndices(u_int32_t nVars)
{
    std::vector<u_int32_t> indices;
    auto addIndex = [&](u_int32_t idx) {
        if (idx < nVars) {
            indices.push_back(idx);
        }
    };
    auto addRange = [&](u_int32_t start, u_int32_t end) {
        if (start >= nVars) {
            return;
        }
        const u_int32_t cappedEnd = std::min(end, static_cast<u_int32_t>(nVars - 1));
        for (u_int32_t idx = start; idx <= cappedEnd; ++idx) {
            indices.push_back(idx);
        }
    };

    addIndex(1);
    addIndex(13);
    addIndex(127);
    addIndex(165);
    addRange(5027, 5266);
    addRange(5271, 5278);
    addRange(5516, 6226);

    return indices;
}

std::vector<u_int32_t> BuildCommonIndices(u_int32_t nVars, const std::vector<u_int32_t> &modifiedIndices)
{
    std::vector<bool> isModified(nVars, false);
    for (u_int32_t idx : modifiedIndices) {
        if (idx < nVars) {
            isModified[idx] = true;
        }
    }

    std::vector<u_int32_t> indices;
    indices.reserve(nVars - std::min<u_int32_t>(nVars, modifiedIndices.size()));
    for (u_int32_t idx = 0; idx < nVars; ++idx) {
        if (!isModified[idx]) {
            indices.push_back(idx);
        }
    }

    return indices;
}

std::string BuildPublicString(AltBn128::FrElement *wtnsData, uint32_t nPublic)
{
    json jsonPublic;
    AltBn128::FrElement aux;
    for (u_int32_t i = 1; i <= nPublic; i++) {
        AltBn128::Fr.toMontgomery(aux, wtnsData[i]);
        jsonPublic.push_back(AltBn128::Fr.toString(aux));
    }

    return jsonPublic.dump();
}

class PoqBatchProver {
public:
    explicit PoqBatchProver(const std::string &zkeyFilename)
        : zkey(zkeyFilename, "zkey", 1),
          zkeyHeader(ZKeyUtils::loadHeader(&zkey))
    {
        if (!PrimeIsValid(zkeyHeader->rPrime)) {
            throw std::invalid_argument("zkey curve not supported");
        }

        prover = Groth16::makeProver<AltBn128::Engine>(
            zkeyHeader->nVars,
            zkeyHeader->nPublic,
            zkeyHeader->domainSize,
            zkeyHeader->nCoefs,
            zkeyHeader->vk_alpha1,
            zkeyHeader->vk_beta1,
            zkeyHeader->vk_beta2,
            zkeyHeader->vk_delta1,
            zkeyHeader->vk_delta2,
            zkey.getSectionData(4),
            zkey.getSectionData(5),
            zkey.getSectionData(6),
            zkey.getSectionData(7),
            zkey.getSectionData(8),
            zkey.getSectionData(9));

        modifiedIndices = BuildModifiedIndices(zkeyHeader->nVars);
        std::sort(modifiedIndices.begin(), modifiedIndices.end());
        modifiedIndices.erase(std::unique(modifiedIndices.begin(), modifiedIndices.end()), modifiedIndices.end());
        commonIndices = BuildCommonIndices(zkeyHeader->nVars, modifiedIndices);
    }

    void setBaseWitness(BinFileUtils::BinFile &wtnsFile)
    {
        auto wtnsHeader = WtnsUtils::loadHeader(&wtnsFile);
        ValidateWitnessHeader(wtnsHeader.get());

        auto wtnsData = static_cast<AltBn128::FrElement *>(wtnsFile.getSectionData(2));

        prover->computeMSMForIndices(wtnsData, commonIndices, commonA, commonB1, commonB2, commonC);

        hasCommon = true;
    }

    std::unique_ptr<Groth16::Proof<AltBn128::Engine>> proveFullWitness(
        BinFileUtils::BinFile &wtnsFile,
        std::string &publicOutput)
    {
        auto wtnsHeader = WtnsUtils::loadHeader(&wtnsFile);
        ValidateWitnessHeader(wtnsHeader.get());

        auto wtnsData = static_cast<AltBn128::FrElement *>(wtnsFile.getSectionData(2));
        publicOutput = BuildPublicString(wtnsData, zkeyHeader->nPublic);
        return prover->prove(wtnsData);
    }

    std::unique_ptr<Groth16::Proof<AltBn128::Engine>> proveWitness(
        BinFileUtils::BinFile &wtnsFile,
        std::string &publicOutput)
    {
        if (!hasCommon) {
            throw std::runtime_error("Base witness not initialized");
        }

        auto wtnsHeader = WtnsUtils::loadHeader(&wtnsFile);
        ValidateWitnessHeader(wtnsHeader.get());

        auto wtnsData = static_cast<AltBn128::FrElement *>(wtnsFile.getSectionData(2));

        typename AltBn128::Engine::G1Point variableA;
        typename AltBn128::Engine::G1Point variableB1;
        typename AltBn128::Engine::G2Point variableB2;
        typename AltBn128::Engine::G1Point variableC;
        prover->computeMSMForIndices(wtnsData, modifiedIndices, variableA, variableB1, variableB2, variableC);

        AltBn128::Engine &E = AltBn128::Engine::engine;
        typename AltBn128::Engine::G1Point pi_a;
        typename AltBn128::Engine::G1Point pib1;
        typename AltBn128::Engine::G2Point pi_b;
        typename AltBn128::Engine::G1Point pi_c;

        E.g1.add(pi_a, commonA, variableA);
        E.g1.add(pib1, commonB1, variableB1);
        E.g2.add(pi_b, commonB2, variableB2);
        E.g1.add(pi_c, commonC, variableC);

        publicOutput = BuildPublicString(wtnsData, zkeyHeader->nPublic);

        return prover->proveWithPrecomputed(wtnsData, pi_a, pib1, pi_b, pi_c);
    }

private:
    void ValidateWitnessHeader(WtnsUtils::Header *wtnsHeader)
    {
        if (zkeyHeader->nVars != wtnsHeader->nVars) {
            throw std::invalid_argument("Invalid witness length. Circuit: "
                                        + std::to_string(zkeyHeader->nVars)
                                        + ", witness: "
                                        + std::to_string(wtnsHeader->nVars));
        }

        if (!PrimeIsValid(wtnsHeader->prime)) {
            throw std::invalid_argument("different wtns curve");
        }
    }

    BinFileUtils::BinFile zkey;
    std::unique_ptr<ZKeyUtils::Header> zkeyHeader;
    std::unique_ptr<Groth16::Prover<AltBn128::Engine>> prover;
    std::vector<u_int32_t> modifiedIndices;
    std::vector<u_int32_t> commonIndices;
    bool hasCommon = false;

    typename AltBn128::Engine::G1Point commonA;
    typename AltBn128::Engine::G1Point commonB1;
    typename AltBn128::Engine::G2Point commonB2;
    typename AltBn128::Engine::G1Point commonC;
};

} // namespace

int main(int argc, char **argv)
{
    if (argc < 6) {
        std::cerr << "Invalid number of parameters" << std::endl;
        std::cerr << "Usage: poqprover <circuit.zkey> <proof_count> "
                     "<witness_1.wtns> ... <witness_N.wtns> "
                     "<proof_1.json> ... <proof_N.json> "
                     "<public_1.json> ... <public_N.json>" << std::endl;
        return EXIT_FAILURE;
    }

    try {
        const std::string zkeyFilename = argv[1];
        const int proofCount = std::stoi(argv[2]);
        if (proofCount <= 0) {
            throw std::invalid_argument("proof_count must be positive");
        }

        const int expectedArgs = 3 + (proofCount * 3);
        if (argc != expectedArgs) {
            throw std::invalid_argument("Invalid number of parameters for proof_count");
        }

        std::vector<std::string> witnessFiles;
        std::vector<std::string> proofFiles;
        std::vector<std::string> publicFiles;
        witnessFiles.reserve(proofCount);
        proofFiles.reserve(proofCount);
        publicFiles.reserve(proofCount);

        int argIndex = 3;
        for (int i = 0; i < proofCount; ++i) {
            witnessFiles.emplace_back(argv[argIndex++]);
        }
        for (int i = 0; i < proofCount; ++i) {
            proofFiles.emplace_back(argv[argIndex++]);
        }
        for (int i = 0; i < proofCount; ++i) {
            publicFiles.emplace_back(argv[argIndex++]);
        }

        PoqBatchProver batchProver(zkeyFilename);

        for (int i = 0; i < proofCount; ++i) {
            BinFileUtils::BinFile wtnsFile(witnessFiles[i], "wtns", 2);

            if (proofCount > 1 && i == 0) {
                batchProver.setBaseWitness(wtnsFile);
            }

            std::string publicOutput;
            auto proof = (proofCount == 1)
                ? batchProver.proveFullWitness(wtnsFile, publicOutput)
                : batchProver.proveWitness(wtnsFile, publicOutput);
            const std::string proofOutput = proof->toJson().dump();

            std::ofstream proofFile(proofFiles[i]);
            proofFile.write(proofOutput.data(), proofOutput.size());

            std::ofstream publicFile(publicFiles[i]);
            publicFile.write(publicOutput.data(), publicOutput.size());
        }

    } catch (std::exception &e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
