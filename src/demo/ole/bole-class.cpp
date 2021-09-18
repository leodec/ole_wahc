
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

string addr = "localhost";
uint32_t port = 9090;

const ui32 numBlocks = pow(2,7);

const ui32 numProtocols = 10;
const double std_dev = 3.2;


template <typename SchemeType>
void bole_receiver() {

    double start, end;

    //
    // Networking Setup
    //
    cout << "Receiver\n";
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
    cout << "Client key generation: " << (end - start) << " ms\n";

    sendPublicKey(kpSeeded.pkSeeded, chl);
    cout << "Client sent public key\n";
    SecretKey& sk = kpSeeded.sk;

    //
    // Input generation
    //
    using encoding_context_t = typename SchemeType::encoding_context_t;
    using encoding_input_t = typename encoding_context_t::encoding_input_t;

    BOLEReceiverInput<encoding_input_t> receiverInput(numBlocks);
    for (ui32 i = 0; i < numBlocks; i++) 
        receiverInput.x[i] = encoding_context_t::generateRandomInput();
    receiverInput.send(chl);
    

    //
    // BOLE Online
    //

    BOLEReceiverOutput<encoding_input_t> bole_output;

    double totalTime = 0;
    for (ui32 i = 0; i < numProtocols; i++) {
        start = currentDateTime();
        bole_output = BOLEReceiver::online(receiverInput, sk, scheme, chl);
        end = currentDateTime();
        totalTime += (end - start);
    }
    double time = totalTime / numProtocols;
    // cout << "BOLE online time = " << time << " ms\n";
    // cout << "Num OLEs = " << numBlocks*scheme.phim << endl;
    // cout << "per OLE online time = " << 1000*(end - start)/(numBlocks*scheme.phim) << " us\n";

    bole_output.send(chl);

    cout << "BOLE online time = " << time << " ms\n";
    cout << "Num OLEs = " << numBlocks*scheme.phim << endl;
    cout << "per OLE online time = " << 1000*(time)/(numBlocks*scheme.phim) << " us\n";

    chl.close();
    sess.stop();
    ios.stop();
};

template <typename SchemeType>
void bole_sender() {

    // double start, end;

    //
    // Networking Setup
    //
    cout << "Sender\n";
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

    cout << "public key received\n";


    //
    // Data Generation
    //
    using encoding_context_t = typename SchemeType::encoding_context_t;
    using encoding_input_t = typename encoding_context_t::encoding_input_t;

    vector<encoding_input_t> aVecs(numBlocks);
    vector<encoding_input_t> bVecs(numBlocks);
    for (ui32 i = 0; i < numBlocks; i++) {
        aVecs[i] = encoding_context_t::generateRandomInput();
        bVecs[i] = encoding_context_t::generateRandomInput();
    }

    BOLESenderInput<encoding_input_t> senderInput(aVecs, bVecs);

    auto recvIn =  BOLEReceiverInput<encoding_input_t>::receive(chl);

    BOLEReceiverOutput<encoding_input_t> correct = bole_pt<encoding_context_t>(recvIn, senderInput);

    //
    // Run OLE
    //

    cout << "Sender beginning online\n";
    double totalTime = 0;
    double start, end;
    for (ui32 i = 0; i < numProtocols; i++) {
        start = currentDateTime();
        BOLESender::online(senderInput, pk, scheme, chl);
        end = currentDateTime();
        totalTime += end-start;
    }
    double time = totalTime/numProtocols;
    cout << "Server online time = " << time << " ms\n";
    cout << "Server per BOLE time = " << (time)/(scheme.phim*numBlocks) * 1000 << " us\n";

    auto bole_output = BOLEReceiverOutput<encoding_input_t>::receive(chl);

    assert(BOLEReceiverOutput<encoding_input_t>::eq(bole_output, correct));
    cout << "BOLE computed correct result\n";

    chl.close();
    sess.stop();
    ios.stop();
};


template <typename SchemeType>
void launch_ole_batch(int argc, char** argv) {
    if (argc == 1) {
		vector<thread> thrds(2);
		thrds[0] = thread([]() { bole_receiver<SchemeType>(); });
		thrds[1] = thread([]() { bole_sender<SchemeType>(); });

        for (auto& thrd : thrds)
            thrd.join();
	} else if(argc == 2) {
		int role = atoi(argv[1]); // 0: send, 1: recv
		role ? bole_sender<SchemeType>() : bole_receiver<SchemeType>();
	} else if(argc == 3) {
		int role = atoi(argv[1]); // 0: send, 1: recv
		role ? bole_sender<SchemeType>() : bole_receiver<SchemeType>();
	}
    else {
      cout << "this program takes a runtime argument.\n\n"
        << "to run the OLE protocol, run\n\n"
        << "    ole-online [0|1]\n\n"
        << "the optional {0,1} argument specifies in which case the program will\n"
        << "run between two terminals, where each one was set to the opposite value. e.g.\n\n"
        << "    bole-online 0\n\n"
        << "    bole-online 1\n\n"
        << "These programs are fully networked and try to connect at " << addr << ":" << port << ".\n"
        << endl;
    }
};

template <typename ptT, ptT p, ui32 numLimbs>
void run_bole(int argc, char** argv, const bool comm_optimized) {

    constexpr ui32 logn = 13;

    cout << "\n================================================\n";
    cout << "Running BOLE logp="<<log2(p)<<"\n";
    cout << "================================================\n";

    typedef DCRT_Poly_Ring<params<ptT>, logn> PlaintextRing;
    typedef EncodingContext<PlaintextRing, p> encoding_context_t;

    typedef DCRT_Ring<fast_four_limb_reduction_params> IntCryptoRing;
    typedef DCRT_Fast_Four_Limb_Reduction_Params<IntCryptoRing, p> dcrt_params_t;
    typedef BFV_DCRT<encoding_context_t, dcrt_params_t> SchemeType;

    // if (comm_optimized) launch_ole_batch_comm_opt<SchemeType>(argc, argv);
	// else 
    launch_ole_batch<SchemeType>(argc, argv);
};

int main(int argc, char** argv) {
    CHECK_DEBUG_VERBOSE;

    const bool comm_optimized = false;

    // run_bole<ui32, 65537, 4>(argc, argv, comm_optimized);  // 16 bits
    // run_bole<ui32, 1032193, 4>(argc, argv, comm_optimized);  // 20 bits
    // run_bole<ui32, 16760833ULL, 4>(argc, argv, comm_optimized);  // 24 bits
    // run_bole<ui32, 268369921ULL, 4>(argc, argv, comm_optimized);  // 28 bits
    run_bole<ui64, 4294475777ULL, 4>(argc, argv, comm_optimized);  // 32 bits
    // run_bole<ui64, 1099511480321ULL, 4>(argc, argv, comm_optimized);  // 40 bits
}
