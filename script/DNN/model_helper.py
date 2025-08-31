from matplotlib import pyplot as plt

##################################################
# for printing the model structure
def count_parameters(model):
    total = 0
    print(f"{'Layer':<40} {'Param #':>10}")
    print("="*52)
    for name, param in model.named_parameters():
        if param.requires_grad:
            num_params = param.numel()
            total += num_params
            print(f"{name:<40} {num_params:>10}")
    print("="*52)
    print(f"{'Total Trainable Params':<40} {total:>10}")

##################################################
# for plotting the loss change during the training
def plot_loss(train_losses, val_losses, modelName):
    plt.plot(train_losses, label='Train Loss')
    plt.plot(val_losses, label='Val Loss')
    plt.xlabel('Epoch')
    plt.ylabel('Loss')
    plt.legend()
    plt.grid()
    plt.title('Training and Validation Loss')
    plt.savefig(modelName+"_loss.png")
