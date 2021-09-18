
/**
@file ole-pipeline.cpp  --  Benchmark of OLE pipeline stages
*/


#include <stdio.h>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <time.h>

#include <cryptoTools/Common/Defines.h>
#include <cryptoTools/Common/Timer.h>
#include <cryptoTools/Common/Log.h>

#include <cryptoTools/Network/Channel.h>
#include <cryptoTools/Network/Session.h>
#include <cryptoTools/Network/IOService.h>

#include "pke/ole.h"
#include "pke/gazelle-network.h"
#include "utils/debug.h"

using namespace lbcrypto;
using namespace osuCrypto;

std::string addr = "localhost";
uint32_t port = 9090;

const ui32 numBlocks = 1000;
const double std_dev = 3.2;


template <typename SchemeType>
void vole_receiver() {

    double start, end;

    //
    // Networking Setup
    //
    std::cout << "Receiver\n";
    IOService ios;
    Session sess(ios, addr, port, EpMode::Server);
    Channel chl = sess.addChannel();

    //
    // Scheme Setup
    //

    SchemeType scheme(std_dev);

    //
    // Key Setup
    //
    using KeyPairSeeded = typename SchemeType::KeyPairSeeded;
    using SecretKey = typename SchemeType::SecretKey;

    start = currentDateTime();
    KeyPairSeeded kpSeeded = scheme.KeyGenSeeded();
    end = currentDateTime();
    std::cout << "Client key generation: " << (end - start) << " ms\n";

    sendPublicKey(kpSeeded.pkSeeded, chl);
    std::cout << "Client sent public key\n";
    SecretKey& sk = kpSeeded.sk;

    //
    // Input generation
    //
    using encoding_context_t = typename SchemeType::encoding_context_t;
    using IntType = typename encoding_context_t::value_type;

    IntType x = encoding_context_t::generateRandomScalar();
    VOLEReceiverInput<IntType> receiverInput(x);
    chl.send(x);

    //
    // VOLE Online
    //

    start = currentDateTime();
    auto vole_output = VOLEReceiver::online(receiverInput, sk, scheme, chl);
    end = currentDateTime();
    std::cout << "VOLE online time = " << (end - start)/numBlocks << " ms\n";
    std::cout << "per OLE online time = " << 1000*(end - start)/(numBlocks*scheme.phim) << " us\n";

    vole_output.send(chl);
};

template <typename SchemeType>
void vole_sender() {

    double start, end;

    //
    // Networking Setup
    //
    std::cout << "Sender\n";
    IOService ios;
    Session sess(ios, addr, port, EpMode::Client);
    Channel chl = sess.addChannel();

    //
    // Scheme Setup
    //

    const SchemeType scheme(std_dev);

    // Receive public key
    using SeededPublicKey= typename SchemeType::PublicKeySeeded;
    using PublicKey = typename SchemeType::PublicKey;
    SeededPublicKey seededPK;
    receivePublicKey(seededPK, chl);
    PublicKey pk = seededPK.expand();

    std::cout << "public key received\n";


    //
    // Data Generation
    //
    using encoding_context_t = typename SchemeType::encoding_context_t;
    using encoding_input_t = typename encoding_context_t::encoding_input_t;
    using IntType = typename encoding_context_t::value_type;

    auto senderInput = VOLESenderInput<encoding_input_t>::template random<encoding_context_t>(numBlocks);

    IntType x; chl.recv(x);
    VOLEReceiverInput recvIn(x);

    VOLEReceiverOutput<encoding_input_t> correct = vole_pt<encoding_context_t>(recvIn, senderInput);

    //
    // Run OLE
    //

    std::cout << "Sender beginning online\n";
    start = currentDateTime();
    VOLESender::online(senderInput, pk, scheme, chl);
    end = currentDateTime();
    std::cout << "VOLE length = " << scheme.phim*numBlocks << std::endl;
    std::cout << "Server online time = " << (end-start) << " ms\n";
    std::cout << "Server per VOLE time = " << (end-start)/(scheme.phim*numBlocks) * 1000 << " us\n";

    auto vole_output = VOLEReceiverOutput<encoding_input_t>::receive(chl);

    assert(VOLEReceiverOutput<encoding_input_t>::eq(vole_output, correct));
    std::cout << "VOLE computed correct result\n";
};


template <typename SchemeType>
void launch_ole_batch(int argc, char** argv) {
    if (argc == 1) {
		std::vector<std::thread> thrds(2);
		thrds[0] = std::thread([]() { vole_receiver<SchemeType>(); });
		thrds[1] = std::thread([]() { vole_sender<SchemeType>(); });

        for (auto& thrd : thrds)
            thrd.join();
	} else if(argc == 2) {
		int role = atoi(argv[1]); // 0: send, 1: recv
		role ? vole_sender<SchemeType>() : vole_receiver<SchemeType>();
	} else if(argc == 3) {
		int role = atoi(argv[1]); // 0: send, 1: recv
		role ? vole_sender<SchemeType>() : vole_receiver<SchemeType>();
	}
    else {
      std::cout << "this program takes a runtime argument.\n\n"
        << "to run the OLE protocol, run\n\n"
        << "    ole-online [0|1]\n\n"
        << "the optional {0,1} argument specifies in which case the program will\n"
        << "run between two terminals, where each one was set to the opposite value. e.g.\n\n"
        << "    vole-online 0\n\n"
        << "    vole-online 1\n\n"
        << "These programs are fully networked and try to connect at " << addr << ":" << port << ".\n"
        << std::endl;
    }
};

template <typename ptT, ui32 logp, ui32 numLimbs>
void run_vole(int argc, char** argv, const bool comm_optimized) {

    constexpr ui32 logn = 13;
    constexpr ptT p = (1ULL)<<logp;

    std::cout << "\n================================================\n";
    std::cout << "Running VOLE logp="<<logp<<"\n";
    std::cout << "================================================\n";

    typedef DCRT_Poly_Ring<params<ptT>, logn> PlaintextRing;
    typedef NullEncodingContext<PlaintextRing, p> encoding_context_t;
    typedef DCRT_Ring<params<ui64>> IntCryptoRing;
    typedef DCRT_Params<IntCryptoRing, numLimbs, 0, encoding_context_t::modulus> dcrt_params_t;
    typedef BFV_DCRT<encoding_context_t, dcrt_params_t> SchemeType;

    // if (comm_optimized) launch_ole_batch_comm_opt<SchemeType>(argc, argv);
	// else 
    launch_ole_batch<SchemeType>(argc, argv);
};

int main(int argc, char** argv) {
    CHECK_DEBUG_VERBOSE;

    const bool comm_optimized = false;

    // run_vole<ui32, 8, 3>(argc, argv, comm_optimized);
    // run_vole<ui32, 16, 4>(argc, argv, comm_optimized);
    // run_vole<ui32, 20, 4>(argc, argv, comm_optimized);
    // run_vole<ui32, 24, 4>(argc, argv, comm_optimized);
    // run_vole<ui32, 28, 4>(argc, argv, comm_optimized);
    run_vole<ui64, 32, 4>(argc, argv, comm_optimized);
    // run_vole<ui64, 40, 4>(argc, argv, comm_optimized);

    // run_bole<ui32, 65537, 4>(argc, argv, comm_optimized);
    // run_bole<ui32, 1032193, 4>(argc, argv, comm_optimized);
    // run_bole<ui32, 16760833ULL, 4>(argc, argv, comm_optimized);
    // run_bole<ui32, 268369921ULL, 4>(argc, argv, comm_optimized);
    // run_bole<ui64, 4294475777ULL, 4>(argc, argv, comm_optimized);
    // run_bole<ui64, 1099511480321ULL, 4>(argc, argv, comm_optimized);
}
