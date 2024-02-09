

log_output_file = 'test/4ahu_4ahu.log'
with open(log_output_file,'r') as file:
    lines = [line for line in file.readlines()[23:]]
    predicted_affinity_per_comples =[line[12:17] for line in lines]
    CNN_pose_score = [line[24:30] for line in lines]
    CNN_affinity = [line[36:41] for line in lines]

print(CNN_affinity)
print(CNN_pose_score)
print(predicted_affinity_per_comples)