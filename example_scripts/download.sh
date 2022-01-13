#!/bin/bash
mkdir Airways Gastro Skin Oral Uro
java -jar ./camiClient.jar -d https://openstack.cebitec.uni-bielefeld.de:8080/swift/v1/CAMI_Airways ./Airways -p fq.gz
java -jar ./camiClient.jar -d https://openstack.cebitec.uni-bielefeld.de:8080/swift/v1/CAMI_Oral ./Oral -p fq.gz
java -jar ./camiClient.jar -d https://openstack.cebitec.uni-bielefeld.de:8080/swift/v1/CAMI_Gastrointestinal_tract ./Gastro -p fq.gz
java -jar ./camiClient.jar -d https://openstack.cebitec.uni-bielefeld.de:8080/swift/v1/CAMI_Skin ./Skin -p fq.gz
java -jar ./camiClient.jar -d https://openstack.cebitec.uni-bielefeld.de:8080/swift/v1/CAMI_Urogenital_tract ./Uro -p fq.gz
