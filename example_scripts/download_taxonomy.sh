#!/bin/bash
java -jar ./camiClient.jar -d https://openstack.cebitec.uni-bielefeld.de:8080/swift/v1/CAMI_Airways ./Airways -p taxonomic_profile
java -jar ./camiClient.jar -d https://openstack.cebitec.uni-bielefeld.de:8080/swift/v1/CAMI_Oral ./Oral -p taxonomic_profile
java -jar ./camiClient.jar -d https://openstack.cebitec.uni-bielefeld.de:8080/swift/v1/CAMI_Gastrointestinal_tract ./Gastro -p taxonomic_profile
java -jar ./camiClient.jar -d https://openstack.cebitec.uni-bielefeld.de:8080/swift/v1/CAMI_Skin ./Skin -p taxonomic_profile
java -jar ./camiClient.jar -d https://openstack.cebitec.uni-bielefeld.de:8080/swift/v1/CAMI_Urogenital_tract ./Uro -p taxonomic_profile
