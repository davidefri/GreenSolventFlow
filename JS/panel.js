document.addEventListener('DOMContentLoaded', () => {
    const floatingButton = document.getElementById('floatingButton');
    const floatingPanel = document.getElementById('floatingPanel');
    const closePanelButton = document.getElementById('closePanel');

    // Funzione per aprire il pannello
    floatingButton.addEventListener('click', () => {
        floatingPanel.classList.add('active');
        floatingButton.style.display = 'none'; // Nasconde il pulsante quando il pannello si apre
    });

    // Funzione per chiudere il pannello
    closePanelButton.addEventListener('click', () => {
        floatingPanel.classList.remove('active');
        floatingButton.style.display = 'block'; // Fa riapparire il pulsante quando il pannello si chiude
    });

    // Opzionale: chiudi il pannello cliccando fuori da esso
    document.addEventListener('click', (event) => {
        if (!floatingPanel.contains(event.target) && !floatingButton.contains(event.target) && floatingPanel.classList.contains('active')) {
            floatingPanel.classList.remove('active');
            floatingButton.style.display = 'block'; // Fa riapparire il pulsante anche se si clicca fuori
        }
    });
});